# Rank and Select on Degenerate Strings

This repository contains instructions to reproduce the experiments in the paper ``Rank and Select on Degenerate Strings'' (to be announced). 

The experiments are forked from the experiments from ["Subset Wavelet Trees"](https://doi.org/10.4230/LIPIcs.SEA.2023.4) by Alanko et al., which you can find [here](https://github.com/jnalanko/SubsetWT-Experiments/tree/master). In [another paper](https://doi.org/10.1137/1.9781611977714.20) they provide a number of reductions from rank-select on degenerate strings to rank-select on regular strings. The experiments from that paper can be found [here](https://github.com/jnalanko/SBWT_experiments), but identical implementations of the reductions also appear in the repository we forked. 

For the most part we have not modified their files. The exceptions are: 
- Removed `main.cpp`; contained a benchmark for their rank-pair structures which we did not need.
- Removed `covid_dataset_fasta_headers.txt.gz`, which we did not use. 
- Modified `kmer_search.cpp` to extend the $k$-mer benchmark to also include our structures.
- Modified `include/predecessor.hh` to avoid an overflow bug. 
- Added `size_in_bytes()` to each of the files `Subset*.hh`. 
- Modified the Makefile. 
 
In addition to this we added the files:
- `rank_benchmark.cpp` which benchmarks the subset-rank queries directly, as opposed to indirectly through the SBWT index.
- `include/ReductionTheorem1-iii.hh` which implements Theorem1, reduction 3. 
- `include/DenseSparseDecomposition.hh` which implements the DSD
- `include/WrappedWaveletTrees.hh` which defines a wrapper for SDSL wavelet trees to make them compatible with the rest of the data structures. 
- `include/SIMDRank.hh` which is a rank-select structure using SIMD.  

See the paper for details on how we performed our experiments.

We refer to the README of Alanko et al (see below) for instructions on how to build the project and run the tests, but first we note: 

- The E. Coli data takes 5GB, and the ERR5035349 takes about 18 GB. However, 8GB of RAM should be sufficient to run the experiments. You can limit RAM usage when creating the SBWT index (the algorithm uses disk), and the size of the SBWT indices for each of the data sets is approximately 200MB and 1.2GB, respectively. 
- They refer to the `microbenchmark` build target, which no longer exists. Use the target `rank_query_benchmark` to build the rank benchmark. Run it using `./rank_benchmark index.sbwt`. 
- To verify that your indices are the same as those we tested on, copy `ecoli.sbwt` and `metagenome.sbwt` into the `subset-rank-select` directory and run `md5sum --check checksums`. 


Everything that follows is the original README from Alanko.

---


# SubsetWT experiments

This repository contains instructions and code to reproduce the experiments in the paper "Subset Wavelet Trees" for SEA 2023.

# Downloading the data

Downloading the data requires `curl`, and `fasterq-dump` from the [SRA toolkit](https://hpc.nih.gov/apps/sratoolkit.html).

## 3682 E. coli genomes

```
curl -O https://zenodo.org/record/6577997/files/coli3682_dataset.tar.gz
```

## 17,336,887 metagenomic reads
```
fasterq-dump ERR5035349
```

# Building and running the experiments

First, pull the submodules with:

```
git submodule init
git submodule update
```

Then, go the SBWT submodule and build it using the instructions in the submodule. Then, compile the experiments with:

```
cd sdsl-lite
cd build
cmake .. -DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_C_COMPILER=$(which gcc)
make
cd ../..
make microbenchmark kmer_search_benchmark
```

This creates an executables called `microbenchmark` and `kmer_search`. Both of these take as input a plain-matrix SBWT index. The SBWT index can be built by passing a list of filenames to the SBWT program, one filename per line. For the E. coli genomes, if the genomes are in the directory `coli3682_dataset`, we can create the filename list by running `find coli3682_dataset/ -type f > list.txt`. For the metagenome dataset, the input file list should contain just the line `ERR5035349_1.fastq`. Given a list file, the SBWT index can then be built with:

```
mkdir temp
./SBWT/build/bin/sbwt build -i list.txt -o index.sbwt -k 31 --add-reverse-complements --n-threads 4 --ram-gigas 8 --temp-dir temp
```

This will save the index to index.sbwt. This should not take more than a few hours for either of the datasets in the paper.

To run the microbenchmark of the paper on this index, run `./microbenchmark index.sbwt`. To run the k-mer search benchmark, run `./kmer_search index.sbwt queries.fastq`, where queries.fna is the file containing the queries. In case of the metagenomic read set, we queried a file containing the first 25,000 reads of `ERR5035349_1.fastq` (extract with `head -n 100000 ERR5035349_1.fastq > queries.fastq`), and in case of the E. coli genomes, we queried the genome `coli3682_dataset/GCA_000005845.2_ASM584v2.fna`.



