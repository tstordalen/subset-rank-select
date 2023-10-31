#pragma once

// Wrapper for sdsl wavelet trees to make them compatible
// with the other rank structuRankStructuresres
template<typename wt_type_t>
class WaveletTreeWrapper {

    wt_type_t wt;

    public: 
    WaveletTreeWrapper(){}
    WaveletTreeWrapper(std::vector<char> data){
        auto v = sdsl::int_vector<8>(data.size(), 0,0);
        for(size_t i = 0; i < data.size(); i++) v[i] = data[i];
        sdsl::construct_im(wt, v, 0);
    }

    WaveletTreeWrapper& operator =(const WaveletTreeWrapper& other) {
        if(&other != this){
            wt = other.wt;
            return *this;
        } else return *this; // Assignment to self -> do nothing.
    }

    inline int64_t rank(int64_t pos, char c) const {
        return wt.rank(pos,c);
    }

    int64_t size_in_bytes() const{
        return sdsl::size_in_bytes(wt);
    }
};

// Regular wavelet tree
typedef sdsl::wt_blcd<sdsl::bit_vector,
            sdsl::rank_support_v5<>,
            sdsl::select_support_scan<1>,
            sdsl::select_support_scan<0>> plain_wt_t; // No select sppport

// Compressed wavelet tree
typedef sdsl::wt_blcd<sdsl::rrr_vector<>,
            sdsl::rrr_vector<>::rank_1_type,
            sdsl::rrr_vector<>::select_1_type,
            sdsl::rrr_vector<>::select_0_type> rrr_wt_t;


// Wavelet using bitvectors where rank information is interleaved with 
// the bits in the bitvector, as opposed to an external structure.
typedef sdsl::wt_blcd<sdsl::bit_vector_il<>, 
            sdsl::rank_support_il<>, 
            sdsl::select_support_il<1>, 
            sdsl::select_support_il<0>, 
            sdsl::byte_tree<>> il_wt_t;


// Wrappers for the above wavelet trees 
typedef WaveletTreeWrapper<il_wt_t> wt_il_wrapped_t;
typedef WaveletTreeWrapper<rrr_wt_t> wt_rrr_wrapped_t;
typedef WaveletTreeWrapper<plain_wt_t> wt_plain_wrapped_t;

