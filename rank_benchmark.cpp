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
#include "SubsetMatrixRank.hh"


#include "ReductionTheorem1-iii.hh"
#include "DenseAndSparseDecomposition.hh"
#include "WrappedWaveletTrees.hh"
#include "SIMDRank.hh"

#include <random>

struct Query{
    int64_t pos;
    char symbol;
};

template<typename SBWT_variant_t>
void benchmark_rank(std::string sbwt_index_name, uint64_t n_symbols, uint64_t seed, const vector<Query>& queries, ostream &csv_out, const SBWT_variant_t& sbwt_variant, const string& rs_name){
    // rank structure
    auto rs = sbwt_variant.get_subset_rank_structure();

    int64_t sum = 0;
    int64_t t0 = current_time_micros();
    for(const Query& Q : queries){
        sum += rs.rank(Q.pos, Q.symbol);
    }
    int64_t t1 = current_time_micros();
    csv_out << "rank-benchmark" << ", "
            << sbwt_index_name  << ", "
            << queries.size() << ", "
            << seed << ", "
            << sum  << ", "
            << rs_name << ", "
            << (double)(t1-t0)/queries.size() << ", "                 //us / query
            << (double)(rs.size_in_bytes() * 8) / n_symbols << endl;  //bits / char 

    cerr << "Benchmarking: " << rs_name << endl;
    cerr << "Single rank final sum: " << sum << endl;
    cerr << "Single rank ns/query: " << (double)(t1-t0) / queries.size() * 1000 << endl << endl;
}

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

// Generates a hash value of the queries to verify the same test data was used across
// different runs. The hash function is probably not great, but sufficient for this case. 
std::vector<Query> generate_queries(uint64_t n_sets, uint64_t n_queries, uint64_t seed){
    std::vector<Query> queries(n_queries); 

    std::mt19937_64 generator(seed);
    std::uniform_int_distribution<uint64_t>  r_pos(0, n_sets); 
    std::uniform_int_distribution<uint64_t>  r_alphabet(0, 3);
    char map[4] = {'A', 'C', 'G', 'T'};
    
    for(uint64_t i = 0; i < n_queries; i++){
        queries[i].pos = r_pos(generator); 
        queries[i].symbol = map[r_alphabet(generator)];
    }
    return queries; 
}

template <typename T>
std::pair<T, std::string> variant(std::string name, T variant){
    return {variant, name};
}

int main(int argc, char** argv){

    if(argc == 1){
        cerr << "Please give a plain-matrix sbwt file as input" << endl;
        return 1;
    }

    string sbwt_index_file = string(argv[1]);

    sbwt::throwing_ifstream in(sbwt_index_file, ios::binary);
    string sbwt_variant = load_string(in.stream); // read variant type
    if(sbwt_variant != "plain-matrix"){
        cerr << "Error: input is not a plain-matrix SBWT" << endl;
    }

    sbwt::plain_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);
    cerr << "Loaded a plain matrix SBWT with " << sbwt.number_of_subsets() << " subsets" << endl;
    auto sswt = sbwt.get_subset_rank_structure();

    auto seed = 133742;
    auto n_queries = 20000000; 
    auto n_sets = sswt.A_bits.size();
    auto queries = generate_queries(n_sets, n_queries, seed);

    auto n_symbols_in_index = 0;
    for(size_t i = 0; i < sswt.A_bits.size(); i++){
        n_symbols_in_index += sswt.A_bits[i];
        n_symbols_in_index += sswt.C_bits[i];
        n_symbols_in_index += sswt.G_bits[i];
        n_symbols_in_index += sswt.T_bits[i];
    }

    // From 'Rank and Select on Degenerate Strings', Bille et al.
    // Paper: [TBA]
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
    auto       matrix = variant("Matrix", sbwt);
    auto    concat_ef = variant("Concat (ef)", build_variant<sbwt::mef_concat_sbwt_t>(sbwt));
    auto concat_plain = variant("Concat (plain)", build_variant<sbwt::plain_concat_sbwt_t>(sbwt));
    auto     split_ef = variant("Split (ef)", build_variant<sbwt::mef_split_sbwt_t>(sbwt));
    auto    split_rrr = variant("Split (rrr)", build_variant<sbwt::rrr_split_sbwt_t>(sbwt));
    auto  split_plain = variant("Split (plain)", build_variant<sbwt::plain_split_sbwt_t>(sbwt));
    auto   matrix_rrr = variant("Matrix (rrr)", build_variant<sbwt::rrr_matrix_sbwt_t>(sbwt));

    ofstream results("rank_benchmark_results.csv");
    benchmark_rank(sbwt_index_file, n_symbols_in_index, seed, queries, results,      thm1iii.first,      thm1iii.second); 
    benchmark_rank(sbwt_index_file, n_symbols_in_index, seed, queries, results,      dsd_rrr.first,      dsd_rrr.second); 
    benchmark_rank(sbwt_index_file, n_symbols_in_index, seed, queries, results,     dsd_scan.first,     dsd_scan.second); 
    benchmark_rank(sbwt_index_file, n_symbols_in_index, seed, queries, results,         simd.first,         simd.second); 
    
    benchmark_rank(sbwt_index_file, n_symbols_in_index, seed, queries, results,  swt_rrr_gen.first,  swt_rrr_gen.second); 
    benchmark_rank(sbwt_index_file, n_symbols_in_index, seed, queries, results,     swt_scan.first,     swt_scan.second); 
    benchmark_rank(sbwt_index_file, n_symbols_in_index, seed, queries, results,      swt_rrr.first,      swt_rrr.second); 
    benchmark_rank(sbwt_index_file, n_symbols_in_index, seed, queries, results,    swt_split.first,    swt_split.second); 

    benchmark_rank(sbwt_index_file, n_symbols_in_index, seed, queries, results,       matrix.first,       matrix.second); 
    benchmark_rank(sbwt_index_file, n_symbols_in_index, seed, queries, results,    concat_ef.first,    concat_ef.second); 
    benchmark_rank(sbwt_index_file, n_symbols_in_index, seed, queries, results, concat_plain.first, concat_plain.second); 
    benchmark_rank(sbwt_index_file, n_symbols_in_index, seed, queries, results,     split_ef.first,     split_ef.second); 
    benchmark_rank(sbwt_index_file, n_symbols_in_index, seed, queries, results,    split_rrr.first,    split_rrr.second); 
    benchmark_rank(sbwt_index_file, n_symbols_in_index, seed, queries, results,  split_plain.first,  split_plain.second); 
    benchmark_rank(sbwt_index_file, n_symbols_in_index, seed, queries, results,   matrix_rrr.first,   matrix_rrr.second); 
}