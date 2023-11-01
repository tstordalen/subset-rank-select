// From [ABPV]
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include "SBWT.hh"
#include "variants.hh"
#include "throwing_streams.hh"
#include "RRR_generalization.hh"
#include "SplitStructure.hh"
#include "BitMagic.hh"
#include "SDSL_WT.hh"
#include "NewSubsetWT.hh"
#include "SeqIO.hh"

// Structures from [BGS]
#include "DenseAndSparseDecomposition.hh"
#include "ReductionTheorem1-ii.hh"
#include "WrappedWaveletTrees.hh"

template<typename variant_t>
variant_t build_variant(const sbwt::plain_matrix_sbwt_t& input){
    variant_t var(
        input.get_subset_rank_structure().A_bits,
        input.get_subset_rank_structure().C_bits,
        input.get_subset_rank_structure().G_bits,
        input.get_subset_rank_structure().T_bits,
        input.get_streaming_support(),
        input.get_k(),
        input.number_of_kmers(),
        input.get_precalc_k()
    );
    return var;
}

template<typename variant_t> 
void benchmark(const variant_t& index, const string& index_name, const string& dataset_name, const string& query_file, ostream& csv_out){

    sbwt::SeqIO::Reader<> reader(query_file);

    int64_t sum = 0;

    int64_t total_time_micros = 0;
    int64_t n_kmers = 0;

    while(true){
        int64_t len = reader.get_next_read_to_buffer();
        if(len == 0) break;
        int64_t t0 = current_time_micros();
        //vector<int64_t> results = index.streaming_search(reader.read_buf);
        for(int64_t i = 0; i < len-index.get_k()+1; i++){
            sum += index.search(reader.read_buf + i);
            n_kmers++;
        }
        int64_t t1 = current_time_micros();
        total_time_micros += t1 - t0; 
        //for(int64_t x : results) sum += x;
        //n_kmers += results.size();
    }
    cerr << "Variant: " << index_name << endl;
    cerr << "Time: " << (double)(total_time_micros) / n_kmers << " us/kmer" << endl;
    cerr << "Sum: " << sum << endl;

    int64_t space_bits = index.get_subset_rank_structure().size_in_bytes() * 8;
    //space_bits += index.number_of_subsets(); // streaming support one bit per subset
    space_bits += index.get_precalc().size() * sizeof(pair<int64_t, int64_t>) * 8; // Precalc
    
    csv_out << "kmer-search-test, " 
            << dataset_name << ", " 
            << query_file  << ", "
            << index_name << ", " 
            << (double)(total_time_micros) / n_kmers << ", " 
            << (double)space_bits / index.number_of_kmers() << ", " 
            << sum << endl;
    
    cerr << "Space: " << space_bits << "\n";

    cerr << "Total k-mers queried: " << n_kmers << endl << endl;;
}

int main(int argc, char** argv){

    if(argc == 1){
        cerr << "Please give a plain-matrix sbwt file and a fastq query file as input" << endl;
        return 1;
    } 

    string sbwt_index_file = string(argv[1]);
    string query_file = string(argv[2]);

    sbwt::throwing_ifstream in(sbwt_index_file, ios::binary);
    string variant = load_string(in.stream); // read variant type
    if(variant != "plain-matrix"){
        cerr << "Error: input is not a plain-matrix SBWT" << endl;
    }

    sbwt::plain_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);
    cerr << "Loaded a plain matrix SBWT with " << sbwt.number_of_subsets() << " subsets" << endl;



    ofstream results("kmer_search_results.csv", std::ios_base::app);
    
    benchmark(sbwt, "Matrix", sbwt_index_file, query_file, results);
    
     
    auto thm1_ii =  build_variant<sbwt::SBWT<ReductionThm1ii<wt_il_wrapped_t>>>(sbwt);
    benchmark(thm1_ii, "Thm1(iii)", sbwt_index_file, query_file, results);
    
    auto dsd_wt_il = build_variant<sbwt::SBWT<DenseSparseDecomp<wt_il_wrapped_t>>>(sbwt);
    benchmark(dsd_wt_il, "DSD", sbwt_index_file, query_file, results);    

    auto dsd_wt_rrr = build_variant<sbwt::SBWT<DenseSparseDecomp<wt_rrr_wrapped_t>>>(sbwt);
    benchmark(dsd_wt_rrr, "DSD (rrr)", sbwt_index_file, query_file, results);    
    
    auto dsd_bitmagic = build_variant<sbwt::SBWT<DenseSparseDecomp<BitMagic>>>(sbwt);
    benchmark(dsd_bitmagic, "DSD (scan)", sbwt_index_file, query_file, results);    

    auto rrr_generalization = build_variant<sbwt::SBWT<NewSubsetWT<RRR_Generalization, RRR_Generalization>>>(sbwt);
    benchmark(rrr_generalization, "SWT (rrr gen.)", sbwt_index_file, query_file, results);

    auto bitmagic = build_variant<sbwt::SBWT<NewSubsetWT<BitMagic, BitMagic>>>(sbwt);
    benchmark(bitmagic, "SWT (scan)", sbwt_index_file, query_file, results);

    auto rrr_wt = build_variant<sbwt::SBWT<NewSubsetWT<SDSL_WT<rrr_wt_t>, SDSL_WT<rrr_wt_t>>>>(sbwt);
    benchmark(rrr_wt, "SWT (rrr)", sbwt_index_file, query_file, results);
    
    auto split = build_variant<sbwt::SBWT<NewSubsetWT<SplitStructure<4>, SplitStructure<3>>>>(sbwt);
    benchmark(split, "SWT (split)", sbwt_index_file, query_file, results);
    
}
