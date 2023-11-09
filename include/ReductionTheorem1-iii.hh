#pragma once
#include "DenseAndSparseDecomposition.hh"


template<typename rank_structure_t>
class ReductionThm1iii {

public:

    rank_structure_t rank_elems;

    // Sparse bitvector with additional rank support
    sdsl::sd_vector<> empty;
    sdsl::sd_vector<>::rank_1_type rank_empty;

    // Regular bitvector with select support
    sdsl::bit_vector B;
    sdsl::bit_vector::select_1_type sel_B;

    // cmap['A' & 0b111] == 0, cmap['C' & 0b111] == 1, 
    // cmap['G' & 0b111] == 2, cmap['T' & 0b111] == 3,
    const std::array<uint8_t,8> cmap = {5,0,5,1,3,5,5,2};

    ReductionThm1iii(){}
    ReductionThm1iii(const sdsl::bit_vector& A_bits,
                    const sdsl::bit_vector& C_bits,
                    const sdsl::bit_vector& G_bits,
                    const sdsl::bit_vector& T_bits){
    

        auto n = A_bits.size();
        std::vector<char> data; 
        auto empty_tmp = sdsl::bit_vector(n, 0);

        B = sdsl::bit_vector(n+1, 0);
        uint64_t bi = 0;

        uint64_t tot_count = 0;
        uint64_t delta;
        uint64_t set_count = 0;
        for(size_t i = 0; i < n; i++){
            delta = 0;
            if (A_bits[i]) data.push_back(0), delta++;
            if (C_bits[i]) data.push_back(1), delta++;
            if (G_bits[i]) data.push_back(2), delta++;
            if (T_bits[i]) data.push_back(3), delta++; 
            
            if (delta == 0) empty_tmp[i] = 1;
            else B[bi] = 1, bi += delta, set_count++;
            tot_count += delta;
        }
        B[bi++] = 1;

        rank_elems = rank_structure_t(data);
        empty = sdsl::sd_vector<>(empty_tmp);
        sdsl::util::init_support(sel_B, &B);
        sdsl::util::init_support(rank_empty, &empty);        
    }

    ReductionThm1iii(const ReductionThm1iii& other){
        assert(&other != this); 
        operator=(other);
    }

    ReductionThm1iii& operator =(const ReductionThm1iii& other) {
        if(&other != this){

            rank_elems = other.rank_elems;

            B = other.B;
            sel_B = other.sel_B;
            sel_B.set_vector(&B);

            empty = other.empty;
            rank_empty = other.rank_empty;
            rank_empty.set_vector(&empty);
            return *this;
        } else return *this; // Assignment to self -> do nothing.
    }

    // Rank of symbol in half-open interval [0..pos)
    int64_t rank(int64_t pos, char c) const {
        c = cmap[c & 0b111]; 
        pos -= rank_empty(pos);
        auto i = sel_B(pos + 1);
        return rank_elems.rank(i, c);
    }
    
    int64_t size_in_bytes() const{
        return rank_elems.size_in_bytes()
             + sdsl::size_in_bytes(B)
             + sdsl::size_in_bytes(sel_B)
             + sdsl::size_in_bytes(empty)
             + sdsl::size_in_bytes(rank_empty);
    }

    int64_t serialize(ostream& os) const{
        std::cerr << "Serialize has not been implemented for ReductionThm1iii\n";
        exit(1);
    }

    void load(istream& is){
        std::cerr << "Load has not been implemented for ReductionThm1iii\n";
        exit(1);
    }
};
