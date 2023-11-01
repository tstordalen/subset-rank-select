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


#include "ReductionTheorem1-ii.hh"
#include "DenseAndSparseDecomposition.hh"
#include "WrappedWaveletTrees.hh"

#include <random>




// --------------------------- IMPORTANT NOTE!! ----------------------------------
// When compiling the rank benchmark, the macro RUNNING_RANK_BENCHMARK must
// be defined if you want to test the simple SubsetMatrixRank structure.
// Compile with -DRUNNING_RANK_BENCHMARK. (the target 'rank_query_benchmark' in the 
// makefile already does this). 
//
// This is because SubsetMatrixRank expects the alphabet {'A','C','G'.'T'} while this test
// uses the alphabet {0,1,2,3}. We made the minimal change to SubsetMatrixRank
// We made minimal changes that will not affect the performance of SubsetMatrixRank
// during other tests, such as kmer_search

struct Query{
    int64_t pos;
    char symbol;
};

template<typename subset_rank_structure_t>
void benchmark_rank(uint64_t n_sets, uint64_t seed, uint64_t queries_id, const vector<Query>& queries, ostream &csv_out, const subset_rank_structure_t& rs, const string& rs_name){
    int64_t sum = 0;
    int64_t t0 = current_time_micros();
    for(const Query& Q : queries){
        sum += rs.rank(Q.pos, Q.symbol);
    }
    int64_t t1 = current_time_micros();
    csv_out << "rank-benchmark" << ", "
            << seed << ", "
            << queries_id << ", " //use this during data processing to verify the same test was used
            << rs_name << ", "
            << (double)(t1-t0)/queries.size() << ", "          //us/query
            << (double)(rs.size_in_bytes() * 8) / n_sets  //bits/char 
            << sum << endl; // use this during data processing to verify that the results are correct

    cerr << "Benchmarking: " << rs_name << endl;
    cerr << "Single rank final sum: " << sum << endl;
    cerr << "Single rank ns/query: " << (double)(t1-t0) / queries.size() * 1000 << endl << endl;
}

// Build subset rank structure
template<typename rank_structure_t>
rank_structure_t build_srs(const sbwt::plain_matrix_sbwt_t& sbwt){
    return rank_structure_t(
                sbwt.get_subset_rank_structure().A_bits, 
                sbwt.get_subset_rank_structure().C_bits,
                sbwt.get_subset_rank_structure().G_bits,
                sbwt.get_subset_rank_structure().T_bits
    );
}


// Generates a hash value of the queries to verify the same test data was used across
// different runs. The hash function is probably not great, but sufficient for this case. 
std::pair<std::vector<Query>, uint64_t> generate_queries(uint64_t n_sets, uint64_t n_queries, uint64_t seed){
    std::vector<Query> queries(n_queries); 

    std::mt19937_64 generator(seed);
    std::uniform_int_distribution<uint64_t>  r_pos(0, n_sets); 
    std::uniform_int_distribution<uint64_t>  r_alphabet(0, 3);
    
    uint64_t hash = 1;
    // values to compute the hash function
    uint64_t a = 9163009843618600053;
    uint64_t b = 264556097;
    uint64_t c = 796964356;
    uint64_t p = 2147483647;
    uint64_t tmp;
    for(uint64_t i = 0; i < n_queries; i++){
        queries[i].pos = r_pos(generator); 
        queries[i].symbol = r_alphabet(generator);

        tmp = a*(queries[i].pos + queries[i].symbol * c) >> (64 - 8); //multiply-shift (ish) to bytes
        hash = (b * hash + tmp) % p; // hashing for sequences 
    }
    
    return {queries, hash};
}

int main(int argc, char** argv){

    if(argc == 1){
        cerr << "Please give a plain-matrix sbwt file as input" << endl;
        return 1;
    }

    if (argc == 2){
        cerr << "Please give a seed for generating queries" << endl;
        return 1;
    }

    if (argc == 3){
        cerr << "Please provide the number of queries to run" << endl;
        return 1;
    }

    auto seed = std::stoi(string(argv[2]));
    auto n_queries = std::stoi(string(argv[3]));

    string sbwt_index_file = string(argv[1]);

    sbwt::throwing_ifstream in(sbwt_index_file, ios::binary);
    string variant = load_string(in.stream); // read variant type
    if(variant != "plain-matrix"){
        cerr << "Error: input is not a plain-matrix SBWT" << endl;
    }

    #ifndef RUNNING_RANK_BENCHMARK
    std::cerr << endl << endl
              << "WARNING: The compiler flag RUNNING_RANK_BENCHMARK is not defined. If you are\n"
              << "         running the rank benchmark for the SubsetMatrixRank structure, you will\n"
              << "         get incorrect results. See the top of rank_benchmark.cpp for an explanation.\n\n"; 
    #endif


    sbwt::plain_matrix_sbwt_t sbwt;
    sbwt.load(in.stream);
    cerr << "Loaded a plain matrix SBWT with " << sbwt.number_of_subsets() << " subsets" << endl;
    auto sswt = sbwt.get_subset_rank_structure();

    auto n_sets = sswt.A_bits.size();

    auto [queries,hash] = generate_queries(n_sets, n_queries, seed);

    // Structures from  [TBD], Bille, Gørtz, Stordalen. We have to specify that the alphabet is 
    // already {0,1,2,3}. The default is to assume {A,C,G,T} and convert it manually. 

    ofstream results("rank_benchmark_results.csv", std::ios_base::app);

    auto simple       = sbwt.get_subset_rank_structure();
    benchmark_rank(n_sets, seed, hash, queries, results, simple     , "Matrix");
    
    auto thm1_ii      = build_srs<ReductionThm1ii<wt_il_wrapped_t, AlphabetType::Reduced_0123>>(sbwt);
    benchmark_rank(n_sets, seed, hash, queries, results, thm1_ii     , "Thm1(iii)");
    //
    auto dsd_wt_il    = build_srs<DenseSparseDecomp<wt_il_wrapped_t, AlphabetType::Reduced_0123>>(sbwt);
    benchmark_rank(n_sets, seed, hash, queries, results, dsd_wt_il   , "DSD");
    //
    auto dsd_wt_rrr   = build_srs<DenseSparseDecomp<wt_rrr_wrapped_t, AlphabetType::Reduced_0123>>(sbwt);
    benchmark_rank(n_sets, seed, hash, queries, results, dsd_wt_rrr  , "DSD (rrr)");
    //
    auto dsd_bitmagic = build_srs<DenseSparseDecomp<BitMagic, AlphabetType::Reduced_0123>>(sbwt);
    benchmark_rank(n_sets, seed, hash, queries, results, dsd_bitmagic, "DSD (scan)");
    //
    // Structures from "Subset Wavelet Trees", Alanko, Biagi, Puglisi, Vuohtoniemi
    //
    auto rrr_gen  =  build_srs<NewSubsetWT<RRR_Generalization, RRR_Generalization>>(sbwt);
    benchmark_rank(n_sets, seed, hash, queries, results, rrr_gen , "SWT (rrr gen.)");
    //
    auto bitmagic =  build_srs<NewSubsetWT<BitMagic, BitMagic>>(sbwt);
    benchmark_rank(n_sets, seed, hash, queries, results, bitmagic, "SWT (scan)");
    //
    auto split    =  build_srs<NewSubsetWT<SplitStructure<4>, SplitStructure<3>>>(sbwt);
    benchmark_rank(n_sets, seed, hash, queries, results, split   , "SWT (split)");
    //
    auto rrr_wt   =  build_srs<NewSubsetWT<SDSL_WT<rrr_wt_t>, SDSL_WT<rrr_wt_t>>>(sbwt);
    benchmark_rank(n_sets, seed, hash, queries, results, rrr_wt  , "SWT (rrr)");
}
