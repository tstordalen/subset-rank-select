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

#include "DenseAndSparseDecomposition.hh"
#include "ReductionTheorem1-iii.hh"
#include "WrappedWaveletTrees.hh"
#include "SIMDRank.hh"

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
        for(int64_t i = 0; i < len-index.get_k()+1; i++){
            sum += index.search(reader.read_buf + i);
            n_kmers++;
        }
        int64_t t1 = current_time_micros();
        total_time_micros += t1 - t0; 
    }
    cerr << "Variant: " << index_name << endl;
    cerr << "Time: " << (double)(total_time_micros) / n_kmers << " us/kmer" << endl;
    cerr << "Sum: " << sum << endl;

    int64_t space_bits = index.get_subset_rank_structure().size_in_bytes() * 8;
    space_bits += index.get_precalc().size() * sizeof(pair<int64_t, int64_t>) * 8; // Precalc
    
    csv_out << "kmer-search-test, " 
            << dataset_name << ", " 
            << query_file  << ", "
            << index_name << ", " 
	        << sum << ", "
            << (double)(total_time_micros) / n_kmers << ", " 
            << (double)space_bits / index.number_of_kmers() << endl;
    
    cerr << "Space: " << space_bits << "\n";

    cerr << "Total k-mers queried: " << n_kmers << endl << endl;;
}

template <typename T>
std::pair<T, std::string> variant(std::string name, T variant){
    return {variant, name};
}

int main(int argc, char** argv){

    if(argc == 1){
        cerr << "Please give a plain-matrix sbwt file and a fastq query file as input" << endl;
        return 1;
    } 

    string sbwt_index_file = string(argv[1]);
    string query_file = string(argv[2]);

    sbwt::throwing_ifstream in(sbwt_index_file, ios::binary);
    string sbwt_variant = load_string(in.stream); // read variant type
    if(sbwt_variant != "plain-matrix"){
        cerr << "Error: input is not a plain-matrix SBWT" << endl;
    }

    sbwt::plain_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);
    cerr << "Loaded a plain matrix SBWT with " << sbwt.number_of_subsets() << " subsets" << endl;
    
    // From 'Rank and Select on Degenerate Strings', Bille et al.
    // Paper: [TBA]
    auto       matrix = variant("Matrix", sbwt);
    auto      thm1iii = variant("Thm1(iii)", build_variant<sbwt::SBWT<ReductionThm1iii<wt_il_wrapped_t>>>(sbwt));
    auto      dsd_rrr = variant("DSD (rrr)", build_variant<sbwt::SBWT<DenseSparseDecomp<wt_rrr_wrapped_t>>>(sbwt));
    auto     dsd_scan = variant("DSD (scan)", build_variant<sbwt::SBWT<DenseSparseDecomp<BitMagic>>>(sbwt));    
    auto         simd = variant("SIMD", build_variant<sbwt::SBWT<DenseSparseDecomp<SIMDRank<4>>>>(sbwt));

    // From 'Subset Wavelet Trees' by Alanko et al.
    // Paper:  https://doi.org/10.4230/LIPIcs.SEA.2023.4
    auto  swt_rrr_gen = variant("SWT (rrr gen.)", build_variant<sbwt::SBWT<NewSubsetWT<RRR_Generalization, RRR_Generalization>>>(sbwt));
    auto     swt_scan = variant("SWT (scan)", build_variant<sbwt::SBWT<NewSubsetWT<BitMagic, BitMagic>>>(sbwt));
    auto      swt_rrr = variant("SWT (rrr)", build_variant<sbwt::SBWT<NewSubsetWT<SDSL_WT<rrr_wt_t>, SDSL_WT<rrr_wt_t>>>>(sbwt));
    auto    swt_split = variant("SWT (split)", build_variant<sbwt::SBWT<NewSubsetWT<SplitStructure<4>, SplitStructure<3>>>>(sbwt));
    
    // From 'Small Searchable k-Spectra via Subset Rank Queries on the Spectral Burrows-Wheeler Transform', Alanko et al.
    // Paper: https://doi.org/10.1137/1.9781611977714.20
    auto    concat_ef = variant("Concat (ef)", build_variant<sbwt::mef_concat_sbwt_t>(sbwt));
    auto concat_plain = variant("Concat (plain)", build_variant<sbwt::plain_concat_sbwt_t>(sbwt));
    auto     split_ef = variant("Split (ef)", build_variant<sbwt::mef_split_sbwt_t>(sbwt));
    auto    split_rrr = variant("Split (rrr)", build_variant<sbwt::rrr_split_sbwt_t>(sbwt));
    auto  split_plain = variant("Split (plain)", build_variant<sbwt::plain_split_sbwt_t>(sbwt));
    auto   matrix_rrr = variant("Matrix (rrr)", build_variant<sbwt::rrr_matrix_sbwt_t>(sbwt));


    ofstream results("kmer_search_results.csv");
    benchmark(      matrix.first,       matrix.second, sbwt_index_file, query_file, results); 
    benchmark(     thm1iii.first,      thm1iii.second, sbwt_index_file, query_file, results); 
    benchmark(     dsd_rrr.first,      dsd_rrr.second, sbwt_index_file, query_file, results); 
    benchmark(    dsd_scan.first,     dsd_scan.second, sbwt_index_file, query_file, results); 
    benchmark(        simd.first,         simd.second, sbwt_index_file, query_file, results); 
    benchmark( swt_rrr_gen.first,  swt_rrr_gen.second, sbwt_index_file, query_file, results); 
    benchmark(    swt_scan.first,     swt_scan.second, sbwt_index_file, query_file, results); 
    benchmark(     swt_rrr.first,      swt_rrr.second, sbwt_index_file, query_file, results); 
    benchmark(   swt_split.first,    swt_split.second, sbwt_index_file, query_file, results); 
    benchmark(   concat_ef.first,    concat_ef.second, sbwt_index_file, query_file, results); 
    benchmark(concat_plain.first, concat_plain.second, sbwt_index_file, query_file, results); 
    benchmark(    split_ef.first,     split_ef.second, sbwt_index_file, query_file, results); 
    benchmark(   split_rrr.first,    split_rrr.second, sbwt_index_file, query_file, results); 
    benchmark( split_plain.first,  split_plain.second, sbwt_index_file, query_file, results); 
    benchmark(  matrix_rrr.first,   matrix_rrr.second, sbwt_index_file, query_file, results); 
}