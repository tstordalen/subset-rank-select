kmer_search_benchmark:
	g++-10 kmer_search.cpp src/SeqIO.cpp src/globals.cpp ./sdsl-lite/build/lib/libsdsl.a -std=c++20 -I sdsl-lite/include/ -O3 -I include -I ./sdsl-lite/build/external/libdivsufsort/include/ -g -o kmer_search -Wno-deprecated-declarations -march=native -DNDEBUG

rank_query_benchmark:
	g++-10 rank_benchmark.cpp src/SeqIO.cpp src/globals.cpp ./sdsl-lite/build/lib/libsdsl.a -std=c++20 -I sdsl-lite/include/ -O3 -I include -I ./sdsl-lite/build/external/libdivsufsort/include/ -g -o rank_benchmark -Wno-deprecated-declarations -march=native -DNDEBUG -DRUNNING_RANK_BENCHMARK
