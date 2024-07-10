#pragma once
#include <iterator>
#include <stdexcept>

template <std::forward_iterator Iter>
class BitReader {
public:
    BitReader(Iter begin, Iter end) : cur_(begin), end_(end) {
    }

    bool ReadBit() {
        if (cur_ == end_) {
            throw std::runtime_error("BitReader::ReadBit: Attempt to read unexisting bit");
        }

        bool bit = ((*cur_) >> (7 - cur_bit_)) & 1;

        ++cur_bit_;

        if (cur_bit_ == 8) {
            ++cur_;
            cur_bit_ = 0;
        }

        return bit;
    }

    uint32_t ReadBits(uint8_t n) {
        if (n > 32) {
            throw std::invalid_argument("BitReader::ReadBits: Requested too many bits");
        }

        uint32_t res = 0;

        for (uint8_t i = 0; i < n; ++i) {
            if (ReadBit()) {
                res |= (1ul << (n - 1 - i));
            }
        }

        return res;
    }

private:
    uint8_t cur_bit_ = 0;
    Iter cur_;
    Iter end_;
};
