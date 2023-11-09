#pragma once

template<typename rank_structure_t>
class DenseSparseDecomp {

public:
    
    // cmap['A' & 0b111] == 0, cmap['C' & 0b111] == 1, 
    // cmap['G' & 0b111] == 2, cmap['T' & 0b111] == 3
    const std::array<uint8_t,8> cmap = {5,0,5,1,3,5,5,2};
    
    // Main rank structure
    rank_structure_t r_elems;

    // Sparse vectors for the alpahbet {0,1,2,3} (i.e., {'A','C','G','T'}) with
    // rank functionality for the extra characters that don't fit in r_elems
    std::array<sdsl::sd_vector<>, 4> aux;
    std::array<sdsl::sd_vector<>::rank_1_type, 4> r_aux;

    // Sparse vector with rank support to count the number of empty sets
    sdsl::sd_vector<> empty_sets;
    sdsl::sd_vector<>::rank_1_type r_empty_sets;

    DenseSparseDecomp(){}
    DenseSparseDecomp(const sdsl::bit_vector& A_bits,
                    const sdsl::bit_vector& C_bits,
                    const sdsl::bit_vector& G_bits,
                    const sdsl::bit_vector& T_bits){

        
        auto n = A_bits.size();
        
        std::vector<char> data; 

        auto tmp_empty = sdsl::bit_vector(n, 0); 
        char bases[4] = {'A', 'C', 'G', 'T'};
        std::array<sdsl::bit_vector, 4> tmps = {
            sdsl::bit_vector(n, 0),
            sdsl::bit_vector(n, 0),
            sdsl::bit_vector(n, 0),
            sdsl::bit_vector(n, 0) 
        };

        std::vector<uint8_t> present;
        for(size_t i = 0; i < n; i++){
            if(A_bits[i]) present.push_back(0);
            if(C_bits[i]) present.push_back(1);
            if(G_bits[i]) present.push_back(2);
            if(T_bits[i]) present.push_back(3);

            if (present.size() == 0) tmp_empty[i] = 1;
            else { 
                // add one to 'data', the rest to the bit vectors
                auto c_index = present.back();
                present.pop_back();
                data.push_back(c_index); 
                while(present.size() > 0) {
                    tmps[present.back()][i] = 1;
                    present.pop_back();
                }
            }
        }
        r_elems = rank_structure_t(data);

        for(int i = 0; i < 4; i++){
            aux[i] = sdsl::sd_vector<>(tmps[i]);
            sdsl::util::init_support(r_aux[i], &aux[i]);
        }
        empty_sets = sdsl::sd_vector<>(tmp_empty);
	sdsl::util::init_support(r_empty_sets, &empty_sets);
    }

    DenseSparseDecomp(const DenseSparseDecomp& other){
        assert(&other != this); 
        operator=(other);
    }

    DenseSparseDecomp& operator =(const DenseSparseDecomp& other) {
        if(&other != this){
            aux = other.aux;
            r_aux = other.r_aux;
            for(int i = 0; i < 4; i++) r_aux[i].set_vector(&aux[i]);
            empty_sets = other.empty_sets;
            r_empty_sets = other.r_empty_sets;
            r_empty_sets.set_vector(&empty_sets);
            r_elems = other.r_elems;
            return *this;
        } else return *this; // Assignment to self -> do nothing.
    }

    // Rank of symbol in half-open interval [0..pos)
    int64_t rank(int64_t pos, char c) const {
        auto zeroes = r_empty_sets(pos);
        c = cmap[c & 0b111];
        auto res = r_elems.rank(pos-zeroes, c);
        return res + r_aux[c](pos);
    }
    
    int64_t size_in_bytes() const{
        auto size = r_elems.size_in_bytes();
        for(int i = 0; i < 4; i++) {
            size += sdsl::size_in_bytes(aux[i]);
            size += sdsl::size_in_bytes(r_aux[i]);
        }
        size += sizeof(cmap);
        size += sdsl::size_in_bytes(empty_sets);
        size += sdsl::size_in_bytes(r_empty_sets);
        return size;
    }

    int64_t serialize(ostream& os) const{
        std::cerr << "Serialize has not been implemented for DenseSparseDecomp\n";
        exit(1);
    }

    void load(istream& is){
        std::cerr << "Load has not been implemented for DenseSparseDecomp\n";
        exit(1);
    }
};


