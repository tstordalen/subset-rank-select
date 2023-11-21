#pragma once

#include <cstdlib>
#include <immintrin.h>
#include <bitset>
#include <cstdint>
#include <vector>



template <uint64_t vectors_per_block>
class FasterSIMDRank {

public:

    static const uint64_t bits_per_base = 2;
    static const uint64_t bits_per_word = 64;
    static const uint64_t bits_per_vector = 512;
    static const uint64_t words_per_vector = bits_per_vector / bits_per_word;
    static const uint64_t all_ones_64 = 0xFFFFFFFFFFFFFFFF;
    static const uint64_t alphabet_size = 4;

    struct LoHi {
        uint64_t lo[words_per_vector];
        uint64_t hi[words_per_vector];
    };

    // The last two words of both lo and hi in the last LoHi will contain the count up to the block
    static const uint64_t bases_per_block = bits_per_vector * vectors_per_block - (alphabet_size/2) * bits_per_word;
    struct Block {
        LoHi vectors[vectors_per_block];

        void set_count(uint64_t count[4]){
            const uint64_t n_lohi = sizeof(vectors)/sizeof(vectors[0]);
            LoHi *last_lohi = vectors + n_lohi - 1;
            uint64_t *last_word_in_lo = ((uint64_t*)last_lohi) + words_per_vector - 1;
            uint64_t *last_word_in_hi = last_word_in_lo + words_per_vector;

            *(last_word_in_lo-1) = count[0];
            *last_word_in_lo = count[1];
            *(last_word_in_hi-1) = count[2];
            *last_word_in_hi = count[3];
        }

        uint64_t get_count(char c) const{
            const uint64_t offsets[4] = {
                sizeof(Block)/sizeof(uint64_t) - words_per_vector - 2,
                sizeof(Block)/sizeof(uint64_t) - words_per_vector - 1,
                sizeof(Block)/sizeof(uint64_t) - 2,
                sizeof(Block)/sizeof(uint64_t) - 1
            };
            return ((uint64_t*)vectors)[offsets[c]];
        }

        // slow utility function for building the data structure
        inline void set_base(uint64_t index, char c){
            if (index >= bases_per_block || c > 3)
                std::cerr << "Error in set_base in block. Exiting\n", exit(1);
            
            const uint64_t lo_hi_index = index/bits_per_vector;
            const uint64_t word_index = (index % bits_per_vector) / bits_per_word;
            const uint64_t index_in_word = index % bits_per_word;
            const uint64_t lo_bit = c & 0b1;
            const uint64_t hi_bit = (c >> 1) & 0b1;
            vectors[lo_hi_index].lo[word_index] |= lo_bit << index_in_word;
            vectors[lo_hi_index].hi[word_index] |= hi_bit << index_in_word;
        }
    };

    Block *blocks; 
    uint64_t n_blocks;
    uint8_t mask[2*vectors_per_block - 1];
    uint8_t *mask_midpoint;

    static void printmm512(__m512i x){
        uint64_t tmp[8];
        memcpy(tmp, &x, sizeof(x));
        std::cerr << "vec[" << tmp[0];
        for(int i = 1; i < 8; i++){
            std::cerr << ", " << tmp[i];
        }
        std::cerr << "]\n";
    }

    FasterSIMDRank(){}

    // expects a sequence over {0,1,2,3}
    FasterSIMDRank(std::vector<char> input){
        if (sizeof(Block) % sizeof(__m512i) != 0) std::cerr << "THE SIZE IS WRONG!\n";

        for(int i = 0; i < 2*vectors_per_block - 1; i++){
            if (i < vectors_per_block) mask[i] = 0xFF;
            else mask[i] = 0; 
        }
        mask_midpoint  = mask + (vectors_per_block - 1); 

        const uint64_t n_bases = input.size();
        n_blocks = (n_bases + bases_per_block - 1) / bases_per_block;
        blocks = (Block*)std::aligned_alloc(sizeof(__m512i), n_blocks*sizeof(Block));
        if (!blocks) std::cerr << "Aligned allocation failed in FasterSIMDRank. Exiting\n", exit(1);
        memset(blocks, 0, sizeof(Block)*n_blocks);

        uint64_t count_so_far[4] = {0,0,0,0};
        uint64_t count = 0;
        uint64_t block_id = 0;
        blocks[0].set_count(count_so_far);
        for(auto c : input){
            count_so_far[c]++;
            blocks[block_id].set_base(count, c);
            count++;
            
            if (count == bases_per_block){
                block_id++;
                count = 0;
                blocks[block_id].set_count(count_so_far);
            }
        }
    }

    FasterSIMDRank(const FasterSIMDRank& other){
        assert(&other != this); 
        operator=(other);
    }

    FasterSIMDRank& operator =(const FasterSIMDRank& other) {
        if(&other != this){
            n_blocks = other.n_blocks;
            blocks = (Block*)std::aligned_alloc(sizeof(__m512i), n_blocks*sizeof(Block));
            if (!blocks) std::cerr << "Aligned allocation failed while copying FasterSIMDRank. Exiting\n", exit(1);
            memcpy(blocks, other.blocks, sizeof(Block) * n_blocks);
            memcpy(mask,   other.mask,   sizeof(mask));
            mask_midpoint = mask + (vectors_per_block - 1); 
            return *this;
        } else return *this; // Assignment to self -> do nothing.
    }

    ~FasterSIMDRank(){
        if (blocks != nullptr) free(blocks);
    }

    inline uint64_t horizontal_sum(__m512i x) const{
        uint64_t tmp[8] __attribute__((aligned(64)));
        _mm512_store_si512((__m512i*)tmp, x);
        uint64_t count = 0;
        for (int i = 0; i < 8; i++) count += tmp[i];
        return count;
    }   


    template <uint8_t ternary_op_immediate>
    int64_t vectorized_count(const Block *block, const uint64_t word_index) const{
        
        // We set *mask_midpoint and select mask_start such that the mask consists of 
        // 'word_index' 1-bits, followed by only zeroes. From the construction of the object, 
        // there are (vectors_per_block - 1) occurrences of 0xFF to the left of the midpoint, 
        // and the same number of occurrences of 0x00 to the right. 
        constexpr uint8_t bits_per_byte = 8; 
        uint8_t *mask_start = mask_midpoint - (word_index / 8);
        *mask_midpoint = (1 << (word_index % bits_per_byte)) - 1;

        constexpr uint64_t unrolling_size = 8;
        constexpr bool is_valid_unrolling_size = unrolling_size <= vectors_per_block && vectors_per_block % unrolling_size == 0;
        static_assert(is_valid_unrolling_size, "unrolling_size must divide the vectors_per_block\n");

        // Avoid loop carried dependencies, allowing proper unrolling. 
        constexpr uint64_t n_accumulators = 4; 
        __m512i v_count[n_accumulators]; 
        memset(v_count, 0, sizeof(v_count));

        // This may or may not be unrolled by the compiler since vectors_per_block is compile-time known
        for(int i = 0; i < vectors_per_block; i++){ 
            const __m512i lo_load   = _mm512_stream_load_si512((__m512i*)(block->vectors + i) + 0);
            const __m512i hi_load   = _mm512_stream_load_si512((__m512i*)(block->vectors + i) + 1);
            const uint8_t keep_mask = mask_start[i];

            // maskz_ternarylogic_epi64 applies /any/ 3-variable boolean function bitwise to each
            // bit in hi_load, lo_load, lo_load below (lo_load appears twice only because we 
            // need three arguments). It uses keep_mask to determine whether it should keep the result
            // or set the correspondingn 64-bit word to zero. The boolean function is determined as follows. 
            // Consider the truth table over three boolean variables: 
            // A B C|
            // 0 0 0|x
            // 0 0 1|x
            // .....|.
            // 1 1 1|x
            //
            // The ternary_op_immediate describes the results column (the least siginificant bit corresponds to (0,0,0)). 
            // E.g., we use 0b00000011 to count the occurrences of '00' (e.g., A). 
            const __m512i is_match_if_keep = _mm512_maskz_ternarylogic_epi64(keep_mask, hi_load, lo_load, lo_load, ternary_op_immediate);
            
            // Population count per 64-bit word
            const __m512i count_match      = _mm512_popcnt_epi64(is_match_if_keep);

            // Again, avoiding loop carried dependencies
            const uint64_t accumulator_i   = i % n_accumulators;
            v_count[accumulator_i] = _mm512_add_epi64(v_count[accumulator_i], count_match);
        } 

        // This sum /does/ have a loop carried dependency. It can be partially mitigated using a better summing pattern. 
        for(int i = 1; i < n_accumulators; i++){
            v_count[0] = _mm512_add_epi64(v_count[i], v_count[0]);
        }

        return horizontal_sum(v_count[0]);
    }

    int64_t rank(int64_t pos, char c) const{
        if (pos == 0) return 0;
        pos--; // [0,pos] inclusive

        const uint64_t block_index      = pos / bases_per_block;
        const uint64_t pos_in_block     = pos % bases_per_block;
        const Block    *block           = blocks + block_index;
        const uint64_t word_index       = pos_in_block / bits_per_word;
        
        uint64_t count_before_vector_index = 0;
        switch(c){
            case 0: 
            count_before_vector_index = vectorized_count<0b00000011>(block, word_index);
            break;
            case 1:
            count_before_vector_index = vectorized_count<0b00001100>(block, word_index);
            break;
            case 2:
            count_before_vector_index = vectorized_count<0b00110000>(block, word_index);
            break;
            case 3:
            count_before_vector_index = vectorized_count<0b11000000>(block, word_index);
            break;
            default:
            std::cerr << "Unknown character " << (int)c <<". Exiting.\n";
            exit(1);
        }

        const uint64_t count_before_block = block->get_count(c);
        

        // Count in the remaining word
        const uint64_t vector_index     = pos_in_block / bits_per_vector;
        const uint64_t pos_in_vector    = pos_in_block % bits_per_vector;
        const uint64_t word_in_vector   = pos_in_vector / bits_per_word;
        const LoHi *last                = block->vectors + vector_index;
        const uint64_t last_lo          = last->lo[word_in_vector];
        const uint64_t last_hi          = last->hi[word_in_vector];
        
        const uint64_t pos_in_last_word = pos_in_block % bits_per_word;
        uint64_t filter_mask            = ((uint64_t)(1) << pos_in_last_word);
        filter_mask                     = (filter_mask - 1) | filter_mask;

        const uint64_t cmp_lo           =  (c & 0b1)        ? all_ones_64 : 0;
        const uint64_t cmp_hi           =  ((c >> 1) & 0b1) ? all_ones_64 : 0;

        const uint64_t match_lo   = ~(cmp_lo ^ last_lo);
        const uint64_t match_hi   = ~(cmp_hi ^ last_hi);
        const uint64_t match_both = match_lo & match_hi;
        const uint64_t masked     = match_both & filter_mask;
        const uint64_t n_match    = std::popcount(masked);

        return n_match + count_before_block + count_before_vector_index;
    }

    int64_t size_in_bytes() const{
        return n_blocks * sizeof(Block) + sizeof(mask);
    }

    int64_t serialize(ostream& os) const{
        std::cerr << "Serialize has not been implemented for FasterSIMDRank\n";
        exit(1);
    }

    void load(istream& is){
        std::cerr << "Load has not been implemented for FasterSIMDRank\n";
        exit(1);
    }
};
