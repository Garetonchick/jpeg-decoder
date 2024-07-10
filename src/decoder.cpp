#include <decoder.h>
#include <glog/logging.h>
#include <sys/types.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ostream>
#include <stdexcept>
#include <numeric>
#include <cassert>
#include <vector>
#include "fft.h"
#include "huffman.h"
#include "bit_reader.h"
#include "jpeg.h"
#include "parser.h"
// #define ENABLE_LOG

namespace {
template <class T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
    for (const auto& val : v) {
        o << val << " ";
    }

    o << "\n";

    return o;
}

constexpr std::array<std::pair<uint8_t, uint8_t>, 64> GenerateZigzagOrder() {
    std::array<std::pair<uint8_t, uint8_t>, 64> order{};
    size_t step = 0;

    for (size_t diag = 0; diag < 15; ++diag) {
        for (size_t i = 0; i <= std::min(diag, 14 - diag); ++i, ++step) {
            order[step].first = (diag < 8 ? i : (7 - i));
            order[step].second = diag - order[step].first;

            if (diag < 8 ? (!(diag & 1)) : (diag & 1)) {
                std::swap(order[step].first, order[step].second);
            }
        }
    }

    return order;
}

constexpr std::array<std::pair<uint8_t, uint8_t>, 64> kZigzagOrder = GenerateZigzagOrder();

class Decoder {
public:
    Decoder() : dct_calc_input_(64), dct_calc_output_(64) {
        dct_calc_ = std::make_unique<DctCalculator>(8, &dct_calc_input_, &dct_calc_output_);
    }

    Image Decode(const JPEGData& data) {
        if (data.encoding_type != EncodingType::BASELINE) {
            throw std::invalid_argument("Decode: Unsupported encoding type");
        }

        std::vector<std::vector<const Section*>> marker_table(256);

        for (const auto& section : data.sections) {
            marker_table[MarkerToSuffix(section.marker)].push_back(&section);
        }

        auto get_sections = [&marker_table](MarkerType type) {
            return marker_table[MarkerToSuffix(type)];
        };

        auto baseline_sections = get_sections(MarkerType::SOF0);

        if (baseline_sections.size() != 1) {
            throw std::invalid_argument(baseline_sections.empty()
                                            ? "Decode: No baseline DCT section"
                                            : "Decode: Multiple baseline DCT sections");
        }

        SectionProcessSOF0(*baseline_sections.front());

        for (auto section : get_sections(MarkerType::DHT)) {
            SectionProcessDHT(*section);
        }

        for (auto section : get_sections(MarkerType::DQT)) {
            SectionProcessDQT(*section);
        }

        auto sos_sections = get_sections(MarkerType::SOS);

        if (sos_sections.size() != 1) {
            throw std::invalid_argument(sos_sections.empty() ? "Decode: No SOS section"
                                                             : "Decode: Multiple SOS sections");
        }

        SectionProcessSOS(*sos_sections.front());

        auto com_sections = get_sections(MarkerType::COM);

        if (!com_sections.empty()) {
            SectionProcessCOM(*com_sections.front());
        }

        Image rval = {};
        std::swap(result_, rval);

        return rval;
    }

private:
    struct ChannelInfo {
        uint8_t id;
        uint8_t dqt_id;
        uint8_t dht_dc_id;
        uint8_t dht_ac_id;
        uint8_t horizontal_thinning;
        uint8_t vertical_thinning;
    };

    using DQTable = std::vector<std::vector<uint16_t>>;
    using ColorTable = std::vector<std::vector<int16_t>>;

    void SectionProcessCOM(const Section& section) {
        std::string comment;
        comment.reserve(section.length);

        for (uint8_t byte : section.data) {
            comment.push_back(*reinterpret_cast<char*>(&byte));
        }

        result_.SetComment(comment);
    }

    void SectionProcessDQT(const Section& section) {
        // DLOG(INFO) << "\n Zigzag Order: \n";

        // for(auto [x, y] : kZigzagOrder) {
        //     DLOG(INFO) << static_cast<uint32_t>(x) << " " << static_cast<uint32_t>(y) << "\n";
        // }

        size_t offset = 0;

        while (offset < section.data.size()) {
            uint8_t first_byte = section.data[offset++];
            uint8_t value_len = ((first_byte & 0xf0) >> 4) + 1;
            uint8_t id = (first_byte & 0x0f);

            if (value_len != 1 && value_len != 2) {
                throw std::invalid_argument("Decoder: Unsupported value length in DQT section");
            }

            if (section.data.size() < offset + 64ull * value_len) {
                throw std::invalid_argument("Decoder: DQT section has wrong size");
            }

            DQTable table(8, std::vector<uint16_t>(8));

            if (value_len == 1) {
                for (size_t i = 0; i < 64; ++i) {
                    auto [x, y] = kZigzagOrder[i];
                    table[x][y] = section.data[i + offset];
                }
            } else {
                for (size_t i = 0; i < 64; ++i) {
                    auto [x, y] = kZigzagOrder[i];
                    uint16_t first_half = section.data[2 * i + offset];
                    table[x][y] = (first_half << 8) + section.data[2 * i + 1 + offset];
                }
            }

            offset += 64 * value_len;

            if (dqtables_.size() <= id) {
                dqtables_.resize(id + 1);
            }

            dqtables_[id] = table;
        }
    }

    void SectionProcessDHT(const Section& section) {
        size_t offset = 0;

        while (offset < section.data.size()) {
            uint8_t first_byte = section.data[offset++];
            uint8_t coef_class = ((first_byte & 0xf0) >> 4);
            uint8_t id = (first_byte & 0x0f);
            std::vector<uint8_t> code_lengths(16);
            std::vector<uint8_t> values;

            if (coef_class > 1) {
                throw std::invalid_argument("Decoder: Invalid DHT section class");
            }

            if (section.data.size() < offset + 16) {
                throw std::invalid_argument("Decoder: Invalid DHT section size");
            }

            size_t values_cnt = 0;

            for (size_t i = 0; i < code_lengths.size(); ++i) {
                code_lengths[i] = section.data[i + offset];
                values_cnt += code_lengths[i];
            }

            offset += 16;

            if (section.data.size() < offset + values_cnt) {
                throw std::invalid_argument("Decoder: Invalid DHT section size");
            }

            values.reserve(values_cnt);

            for (size_t i = 0; i < values_cnt; ++i) {
                values.push_back(section.data[offset + i]);
            }

            offset += values_cnt;

            std::vector<HuffmanTree>* trees = &(coef_class ? dhtrees_ac_ : dhtrees_dc_);

            if (trees->size() <= id) {
                trees->resize(id + 1);
            }

            try {
                (*trees)[id].Build(code_lengths, values);
            } catch (...) {
                throw std::invalid_argument("Decoder: Failed to build huffman tree from DHT table");
            }
        }
    }

    void SectionProcessSOF0(const Section& section) {
        if (section.data.size() < 6) {
            throw std::invalid_argument("Decoder: Invalid SOF0 section size");
        }

        precision_ = section.data.front();
        height_ = CombineBytes(section.data[1], section.data[2]);
        width_ = CombineBytes(section.data[3], section.data[4]);
        result_.SetSize(width_, height_);
        channels_cnt_ = section.data[5];

        if (std::min(width_, height_) == 0) {
            throw std::invalid_argument("Decoder: Zero sized image");
        }

        if (channels_cnt_ != 3 && channels_cnt_ != 1) {
            throw std::invalid_argument("Decoder: Unsupported channels count");
        }

        if (section.data.size() != channels_cnt_ * 3 + 6) {
            throw std::invalid_argument("Decoder: Invalid SOF0 section size");
        }

        uint8_t mxh = 0;
        uint8_t mxv = 0;

        for (uint8_t i = 0; i < channels_cnt_; ++i) {
            ChannelInfo info;
            info.id = section.data[6 + i * 3];
            uint8_t thinning_byte = section.data[6 + i * 3 + 1];
            info.horizontal_thinning = ((thinning_byte & 0xf0) >> 4);
            info.vertical_thinning = (thinning_byte & 0x0f);
            info.dqt_id = section.data[6 + i * 3 + 2];

            if (!info.horizontal_thinning || !info.vertical_thinning) {
                throw std::invalid_argument("Decoder: Invalid thinning");
            }

            mxh = std::max(mxh, info.horizontal_thinning);
            mxv = std::max(mxv, info.vertical_thinning);

            channels_.emplace_back(std::move(info));
        }

        // for(auto& info : channels_) {
        //     info.horizontal_thinning = mxh / info.horizontal_thinning;
        //     info.vertical_thinning = mxv / info.vertical_thinning;
        // }
    }

    void SectionProcessPreambleSOS(const Section& section) {
        uint8_t channels_cnt = section.data.front();

        if (channels_cnt != channels_cnt_) {
            throw std::invalid_argument("Decode: Channels count mismatch in different sections");
        }

        if (section.data.size() < 4 + channels_cnt * 2 || section.length != 4 + channels_cnt * 2) {
            throw std::invalid_argument("Decode: SOS preamble has wrong size");
        }

        if (section.data[1 + channels_cnt * 2] != 0x00 ||
            section.data[2 + channels_cnt * 2] != 0x3f ||
            section.data[3 + channels_cnt * 2] != 0x00) {
            throw std::invalid_argument("Decode: SOS preamble is invalid for baseline");
        }

        for (uint8_t i = 0; i < channels_cnt; ++i) {
            uint8_t channel_id = section.data[1 + i * 2];
            uint8_t second_byte = section.data[2 + i * 2];
            uint8_t dht_dc_id = ((second_byte & 0xf0) >> 4);
            uint8_t dht_ac_id = (second_byte & 0x0f);

            for (auto& info : channels_) {
                if (info.id == channel_id) {
                    info.dht_ac_id = dht_ac_id;
                    info.dht_dc_id = dht_dc_id;

                    if (dht_ac_id >= dhtrees_ac_.size() || dht_dc_id >= dhtrees_dc_.size()) {
                        throw std::invalid_argument("Decode: Invalid huffman tree id");
                    }

                    break;
                }
            }
        }
    }

    int HuffmanRead(BitReader<std::vector<uint8_t>::const_iterator>& bit_reader,
                    HuffmanTree& tree) {
        int val = 0;

        while (!tree.Move(bit_reader.ReadBit(), val)) {
        }

        return val;
    }

    int ReadValue(BitReader<std::vector<uint8_t>::const_iterator>& bit_reader, uint8_t n_bits) {
        if (!n_bits) {
            return 0;
        }

        int val = bit_reader.ReadBits(n_bits);

        if (!((val >> (n_bits - 1)) & 1)) {
            val = val - (1 << n_bits) + 1;
        }

        return val;
    }

    int ReadDC(BitReader<std::vector<uint8_t>::const_iterator>& bit_reader, HuffmanTree& tree) {
        int n_bits = HuffmanRead(bit_reader, tree);

        if (!n_bits) {
            return 0;
        }

        return ReadValue(bit_reader, n_bits);
    }

    std::pair<uint8_t, int> ReadAC(BitReader<std::vector<uint8_t>::const_iterator>& bit_reader,
                                   HuffmanTree& tree) {
        int n_zeros_and_n_bits = HuffmanRead(bit_reader, tree);
        assert(n_zeros_and_n_bits < 256 && n_zeros_and_n_bits >= 0);
        uint8_t n_zeros = ((n_zeros_and_n_bits & 0xf0) >> 4);
        uint8_t n_bits = (n_zeros_and_n_bits & 0x0f);
        assert(n_zeros <= 64);

        if (!n_zeros_and_n_bits) {
            return {64, 0};
        }

        return {n_zeros, ReadValue(bit_reader, n_bits)};
    }

    void InverseDCT(std::vector<int16_t>* block) {
        for (uint8_t i = 0; i < block->size(); ++i) {
            dct_calc_input_[i] = (*block)[i];
        }

        dct_calc_->Inverse();

        for (uint8_t i = 0; i < block->size(); ++i) {
            (*block)[i] = std::round(dct_calc_output_[i]);
        }
    }

    std::vector<int16_t> DecodeBlock(BitReader<std::vector<uint8_t>::const_iterator>& bit_reader,
                                     HuffmanTree& dht_dc, HuffmanTree& dht_ac, DQTable& dqt,
                                     int16_t& prev_dc) {
        std::vector<int16_t> block(64, 0);

        block[0] = ReadDC(bit_reader, dht_dc) + prev_dc;
        prev_dc = block[0];

        for (size_t i = 1; i < block.size();) {
            auto [n_zeros, val] = ReadAC(bit_reader, dht_ac);
            i += n_zeros;

            if (n_zeros == 64) {
                break;
            }

            auto [x, y] = kZigzagOrder.at(i);
            block[x * 8 + y] = val;
            ++i;
        }

        assert(dqt.size() == 8);

        for (uint8_t i = 0; i < dqt.size(); ++i) {
            assert(dqt[0].size() == 8);

            for (uint8_t j = 0; j < dqt[i].size(); ++j) {
                block[i * 8 + j] *= dqt[i][j];
            }
        }

        InverseDCT(&block);

        for (auto& val : block) {
            val = std::clamp(val + 128, 0, 255);
        }

        return block;
    }

    std::vector<ColorTable> DecodeMCU(BitReader<std::vector<uint8_t>::const_iterator>& bit_reader) {
        std::vector<ColorTable> mcu;
        mcu.reserve(channels_.size());

        size_t channel_idx = 0;

        for (const auto& channel : channels_) {
            mcu.emplace_back(channel.vertical_thinning * 8,
                             std::vector<int16_t>(channel.horizontal_thinning * 8));

            for (size_t i = 0; i < channel.vertical_thinning; ++i) {
                for (size_t j = 0; j < channel.horizontal_thinning; ++j) {
                    auto block =
                        DecodeBlock(bit_reader, dhtrees_dc_.at(channel.dht_dc_id),
                                    dhtrees_ac_.at(channel.dht_ac_id), dqtables_.at(channel.dqt_id),
                                    prev_dcs_.at(channel_idx));

                    for (size_t x = 0; x < 8; ++x) {
                        for (size_t y = 0; y < 8; ++y) {
                            mcu.back()[i * 8 + x][j * 8 + y] = block[x * 8 + y];
                        }
                    }
                }
            }

            ++channel_idx;
        }

#ifdef ENABLE_LOG
        DLOG(INFO) << "Mcu is: "
                   << "\n"
                   << mcu << "\n";
#endif

        return mcu;
    }

    std::vector<ColorTable> GetRawChannels(
        BitReader<std::vector<uint8_t>::const_iterator>& bit_reader) {
        prev_dcs_.resize(channels_.size(), 0);
        size_t mxv_thinning =
            std::max_element(channels_.begin(), channels_.end(),
                             [](const auto& info1, const auto& info2) {
                                 return info1.vertical_thinning < info2.vertical_thinning;
                             })
                ->vertical_thinning;
        size_t mxh_thinning =
            std::max_element(channels_.begin(), channels_.end(),
                             [](const auto& info1, const auto& info2) {
                                 return info1.horizontal_thinning < info2.horizontal_thinning;
                             })
                ->horizontal_thinning;
        size_t mcu_width = mxh_thinning * 8;
        size_t mcu_height = mxv_thinning * 8;

        size_t mcu_hor_cnt = (width_ + mcu_width - 1) / mcu_width;
        size_t mcu_vert_cnt = (height_ + mcu_height - 1) / mcu_height;

        std::vector<ColorTable> ctables(channels_.size(),
                                        ColorTable(height_, std::vector<int16_t>(width_)));

        // DLOG(INFO) << "Number of MCU: " << (mcu_hor_cnt * mcu_vert_cnt) << "\n";

        // try {

        for (size_t i = 0; i < mcu_vert_cnt; ++i) {
            for (size_t j = 0; j < mcu_hor_cnt; ++j) {
                auto mcu = DecodeMCU(bit_reader);

                for (size_t x = i * mcu_height; x < std::min<size_t>((i + 1) * mcu_height, height_);
                     ++x) {
                    for (size_t y = j * mcu_width;
                         y < std::min<size_t>((j + 1) * mcu_width, width_); ++y) {
                        for (size_t cidx = 0; cidx < mcu.size(); ++cidx) {
                            size_t mcu_x = (x - i * mcu_height) /
                                           (mxv_thinning / channels_[cidx].vertical_thinning);
                            size_t mcu_y = (y - j * mcu_width) /
                                           (mxh_thinning / channels_[cidx].horizontal_thinning);
                            ctables[cidx][x][y] = mcu[cidx][mcu_x][mcu_y];
                        }
                    }
                }
            }
        }
        // } catch(...) {
        //     DLOG(INFO) << "Caught exception in GetYCbCr\n";
        // }

#ifdef ENABLE_LOG
        DLOG(INFO) << "Decoded image in YCbCr is:\n" << ctables << "\n";
#endif

        return ctables;
    }

    void YCbCrToRGB(int16_t y, int16_t cb, int16_t cr, int16_t* r, int16_t* g, int16_t* b) {
        int16_t tr = std::round(1.0 * y + 1.402 * (cr - 128));
        int16_t tg = std::round(1.0 * y - 0.34414 * (cb - 128) - 0.71414 * (cr - 128));
        int16_t tb = std::round(1.0 * y + 1.772 * (cb - 128));
        // int16_t tr = std::round(1.0 * y + 1.402 * (cr - 128));
        // int16_t tg = std::round(1.0 * y - 0.3 * (cb - 128) - 0.7 * (cr - 128));
        // int16_t tb = std::round(1.0 * y + 1.772 * (cb - 128));
        *r = std::clamp<int16_t>(tr, 0, 255);
        *g = std::clamp<int16_t>(tg, 0, 255);
        *b = std::clamp<int16_t>(tb, 0, 255);
    }

    void ConvertYCbCrToRGB(std::vector<ColorTable>* color_tables) {
        for (size_t i = 0; i < height_; ++i) {
            for (size_t j = 0; j < width_; ++j) {
                YCbCrToRGB((*color_tables)[0][i][j], (*color_tables)[1][i][j],
                           (*color_tables)[2][i][j], &(*color_tables)[0][i][j],
                           &(*color_tables)[1][i][j], &(*color_tables)[2][i][j]);
            }
        }
    }

    std::vector<ColorTable> ConvertGrayscaleToRGB(const std::vector<ColorTable>& color_tables) {
        std::vector<ColorTable> res(3, color_tables[0]);
        for (size_t i = 0; i < height_; ++i) {
            for (size_t j = 0; j < width_; ++j) {
                res[0][i][j] = color_tables[0][i][j];
                res[1][i][j] = color_tables[0][i][j];
                res[2][i][j] = color_tables[0][i][j];
            }
        }

        return res;
    }

    void SectionProcessSOS(const Section& section) {
        SectionProcessPreambleSOS(section);
        BitReader bit_reader(section.data.begin() + section.length, section.data.end());
        auto color_tables = GetRawChannels(bit_reader);

        if (channels_cnt_ == 3) {
            ConvertYCbCrToRGB(&color_tables);
        } else if (channels_cnt_ == 1) {
            color_tables = ConvertGrayscaleToRGB(color_tables);
        }

#ifdef ENABLE_LOG
        DLOG(INFO) << "Decoded image in RGB is:\n" << color_tables << "\n";
#endif

        for (size_t i = 0; i < height_; ++i) {
            for (size_t j = 0; j < width_; ++j) {
                result_.SetPixel(i, j,
                                 RGB{.r = color_tables[0][i][j],
                                     .g = color_tables[1][i][j],
                                     .b = color_tables[2][i][j]});
            }
        }
    }

    uint16_t CombineBytes(uint8_t b1, uint8_t b2) {
        return (static_cast<uint16_t>(b1) << 8) + static_cast<uint16_t>(b2);
    }

private:
    Image result_;
    std::vector<HuffmanTree> dhtrees_dc_;
    std::vector<HuffmanTree> dhtrees_ac_;
    std::vector<DQTable> dqtables_;
    std::vector<ChannelInfo> channels_;
    std::vector<double> dct_calc_input_;
    std::vector<double> dct_calc_output_;
    std::unique_ptr<DctCalculator> dct_calc_;
    std::vector<int16_t> prev_dcs_;
    uint8_t precision_;
    uint8_t channels_cnt_;
    uint16_t height_;
    uint16_t width_;
};
}  // namespace

Image Decode(std::istream& input) {
    JPEGData data = Parse(input);
    Decoder decoder;

    return decoder.Decode(data);
}
