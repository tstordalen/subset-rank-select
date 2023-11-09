kmer_search_benchmark:
	g++-10 -march=native -mavx -mavx512f -mavx512vpopcntdq kmer_search.cpp src/SeqIO.cpp src/globals.cpp ./sdsl-lite/build/lib/libsdsl.a -std=c++20 -I sdsl-lite/include/ -O3 -I include -I ./sdsl-lite/build/external/libdivsufsort/include/ -g -o kmer_search -Wno-deprecated-declarations -DNDEBUG 

rank_query_benchmark:
	g++-10 -march=native -mavx -mavx512f -mavx512vpopcntdq rank_benchmark.cpp src/SeqIO.cpp src/globals.cpp ./sdsl-lite/build/lib/libsdsl.a -std=c++20 -I sdsl-lite/include/ -O3 -I include -I ./sdsl-lite/build/external/libdivsufsort/include/ -g -o rank_benchmark -Wno-deprecated-declarations -DNDEBUG
