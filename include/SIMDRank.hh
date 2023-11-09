#pragma once

#include <cstdlib>
#include <immintrin.h>

template <uint64_t block_size>
class SIMDRank {

public:
    
    const uint64_t bits_per_base = 2;
    const uint64_t bits_per_word = 64;
    const uint64_t bits_per_vector = 512;
    const uint64_t words_per_vector = bits_per_vector / bits_per_word;
    const uint64_t all_ones_64 = 0xFFFFFFFFFFFFFFFF;

    uint64_t *lo_data, *hi_data, *buf;
    uint64_t n_words;

    std::vector<std::array<uint64_t, 3>> checkpoints; 

    // v_comp_constants[0] = all zeroes, v_comp_constants[1] = all ones
    __m512i v_comp_constants[2]; 
    const uint64_t comp_constants[2] = {0, 0xFFFFFFFFFFFFFFFF};
    const uint8_t char_map[4][2] = {
        {0,0},  // for A, pick lo_bit = v_comp_constants[0], hi_bit = v_comp_constants[0]
        {1,0},  // for C, pick lo_bit = v_comp_constants[1], hi_bit = v_comp_constants[0]
        {0,1},  // for G, pick lo_bit = v_comp_constants[0], hi_bit = v_comp_constants[1]
        {1,1}   // for T, pick lo_bit = v_comp_constants[1], hi_bit = v_comp_constants[1]
    };

    __m512i broadcast_u64_to_mm512(uint64_t x){
        __m256i tmp = _mm256_set1_epi64x(x);
        return _mm512_broadcast_i64x4(tmp);
    }

    static void printmm512(__m512i x){
        uint64_t tmp[8];
        memcpy(tmp, &x, sizeof(x));
        std::cerr << "vec[" << tmp[0];
        for(int i = i; i < 8; i++){
            std::cerr << ", " << tmp[i];
        }
        std::cerr << "]\n";
    }

    SIMDRank(){}

    // expects a sequence over {0,1,2,3}
    SIMDRank(std::vector<char> input){
        
        v_comp_constants[0] = broadcast_u64_to_mm512(0);
        v_comp_constants[1] = broadcast_u64_to_mm512(0xFFFFFFFFFFFFFFFF); 

        const uint64_t n_bases = input.size();
        n_words = (n_bases + bits_per_word - 1) / bits_per_word;

        // align data for fast SIMD loads
        lo_data = (uint64_t*)std::aligned_alloc(sizeof(__m512i), n_words * bits_per_word / 8);
        hi_data = (uint64_t*)std::aligned_alloc(sizeof(__m512i), n_words * bits_per_word / 8);
        buf     = (uint64_t*)std::aligned_alloc(sizeof(__m512i), sizeof(__m512i));

        if (!lo_data || !hi_data || !buf) std::cerr << "Aligned allocation failed in SIMDRank. Exiting.\n", exit(1);

        // Take the input stream, separate each character into lo and hi bit
        // and write it to lo_data and hi_data
        uint64_t hi_current = 0;
        uint64_t lo_current = 0;
        uint64_t n_current = 0;
        uint64_t lo_tmp, hi_tmp;
        uint64_t data_i = 0;
        // assumes e in {0,1,2,3}
        for (char e : input){
            lo_tmp = e & 0b1;
            hi_tmp = (e >> 1) & 0b1;
            hi_current |= (hi_tmp << n_current);
            lo_current |= (lo_tmp << n_current);
            n_current++;
            if (n_current == bits_per_word){
                hi_data[data_i] = hi_current;
                lo_data[data_i] = lo_current;
                data_i++;
                hi_current = 0;
                lo_current = 0;
                n_current = 0;
            }
        }
        hi_data[data_i] = hi_current;
        lo_data[data_i] = lo_current;
        data_i++;

        // Compute checkpoints for every 'block_size'th vector
        std::array<uint64_t, 3> count_so_far = {0,0,0};
        checkpoints.push_back(count_so_far);
        uint64_t count = 0;
        for(uint64_t i = 0; i < n_bases; i++){
            if (input[i] < 3) count_so_far[input[i]]++;
            
            count++;

            if(count == block_size * bits_per_vector){
                checkpoints.push_back(count_so_far);
                count = 0;
            }
        }
    }

    SIMDRank(const SIMDRank& other){
        assert(&other != this); 
        operator=(other);
    }

    SIMDRank& operator =(const SIMDRank& other) {
        if(&other != this){
            lo_data = other.lo_data;
            hi_data = other.hi_data;
            buf = other.buf;
            v_comp_constants[0] = other.v_comp_constants[0];
            v_comp_constants[1] = other.v_comp_constants[1];
            checkpoints = other.checkpoints;
            n_words = other.n_words;
            return *this;
        } else return *this; // Assignment to self -> do nothing.
    }

    // Rank of symbol in half-open interval [0..pos)
    int64_t rank(int64_t pos, char c) const {
        if (pos == 0) return 0;
        pos--; //[0,pos] inclusive
            
        uint64_t checkpoint_index = pos / (block_size * bits_per_vector);
        uint64_t count = 0;

        if (c < 3) count = checkpoints[checkpoint_index][c];
        else{
            const uint64_t n_elements_before_checkpoint = checkpoint_index * (block_size * bits_per_vector);
            count =  n_elements_before_checkpoint
                       - checkpoints[checkpoint_index][0]
                       - checkpoints[checkpoint_index][1]
                       - checkpoints[checkpoint_index][2];
        }
        

        // set all fields to zero
        __m512i v_count = v_comp_constants[0];

        // fill v_lo_cmp with the low bit of c, and v_hi_cmp with the high bit of c
        const __m512i v_lo_cmp = v_comp_constants[char_map[c][0]];
        const __m512i v_hi_cmp = v_comp_constants[char_map[c][1]];

        uint64_t simd_current   = checkpoint_index * block_size * bits_per_vector / bits_per_word;
        const uint64_t simd_end = pos / bits_per_vector;
        while (simd_current < simd_end){
            const __m512i lo_load = _mm512_stream_load_si512((__m512i*)(lo_data + simd_current));
            const __m512i hi_load = _mm512_stream_load_si512((__m512i*)(hi_data + simd_current));
            simd_current += words_per_vector;

            const __m512i lo_xor = _mm512_xor_epi64(lo_load, v_lo_cmp);
            const __m512i hi_xor = _mm512_xor_epi64(hi_load, v_hi_cmp);

            const __m512i lo_neg = _mm512_ternarylogic_epi32(lo_xor, lo_xor, lo_xor, 0b01010101); //bitwise negation             
            const __m512i hi_neg = _mm512_ternarylogic_epi32(hi_xor, hi_xor, hi_xor, 0b01010101); //bitwise negation
            const __m512i and_res = _mm512_and_epi64(lo_neg, hi_neg);    

            // accumulate in v_count
            const __m512i cnt = _mm512_popcnt_epi64(and_res);
            v_count = _mm512_add_epi64(v_count, cnt);
        }
        
        // sum over v_count
        _mm512_store_si512((__m512i*)buf, v_count);
        for (int64_t i = 0; i < 8; i++) {
            count += buf[i];
        }
        
       
        // choosing either all zeroes or all ones
        const uint64_t lo_cmp = comp_constants[char_map[c][0]];
        const uint64_t hi_cmp = comp_constants[char_map[c][1]];
 
        uint64_t current = simd_current;
        const uint64_t last_pos = pos / bits_per_word;
        uint64_t lo_tmp, hi_tmp, tmp;
        while (current < last_pos){
            lo_tmp = lo_data[current];
            hi_tmp = hi_data[current];
            current++;

            lo_tmp ^= lo_cmp;
            hi_tmp ^= hi_cmp;
            
            lo_tmp = ~lo_tmp;
            hi_tmp = ~hi_tmp;

            tmp = lo_tmp & hi_tmp;
            
            count += std::popcount(tmp);
        }
        
        const uint64_t n_bases_remaining = pos + 1 - current * bits_per_word;
        if (n_bases_remaining > 0){
            lo_tmp = lo_data[current];
            hi_tmp = hi_data[current];
            current++;

            lo_tmp ^= lo_cmp;
            hi_tmp ^= hi_cmp;
            
            lo_tmp = ~lo_tmp;
            hi_tmp = ~hi_tmp;

            tmp = lo_tmp & hi_tmp;

            // filter out matches beyond p
            tmp <<= bits_per_word - n_bases_remaining;

            count += std::popcount(tmp);
        } 

        return count;
    }
    
    int64_t size_in_bytes() const{
        int64_t size = 0;
        size += sizeof(*lo_data) * n_words;
        size += sizeof(*hi_data) * n_words;
        size += sizeof(n_words);
        size += sizeof(__m512i); // buf
        size += sizeof(uint64_t) * 3 * checkpoints.size();
        size += sizeof(v_comp_constants);
        size += sizeof(comp_constants);
        size += sizeof(char_map);
        return size;
    }

    int64_t serialize(ostream& os) const{
        std::cerr << "Serialize has not been implemented for SIMDRank\n";
        exit(1);
    }

    void load(istream& is){
        std::cerr << "Load has not been implemented for SIMDRank\n";
        exit(1);
    }
};
