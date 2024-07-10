#include "parser.h"
#include <cstdint>
#include <iomanip>
#include <exception>
#include <stdexcept>
#include <glog/logging.h>
#include "jpeg.h"
#include "parser.h"
#include "stream_reader.h"

// #define ENABLE_LOG

JPEGData Parse(std::istream& in) {
    StreamReader reader(in);
    std::vector<Section> sections;
    EncodingType encoding_type = EncodingType::BASELINE;
    bool done = false;

    if (!reader.MustEqual2b(MarkerToVal(MarkerType::SOI))) {
        throw std::invalid_argument("Parse: Invalid jpeg structure");
    }

    while (!done) {
        uint16_t marker_val = reader.Expect2b();

        if (!IsMarker(marker_val)) {
            throw std::invalid_argument("Parse: Expected Marker");
        }

        MarkerType marker = ValToMarker(marker_val);

        if (marker == MarkerType::EOI) {
            throw std::invalid_argument("Parse: EOI before SOS marker");
        }

        if (marker == MarkerType::SOI) {
            throw std::invalid_argument("Parse: Multiple start markers");
        }

        uint16_t len = reader.Expect2b();

        if (len <= 2) {
            throw std::invalid_argument("Parse: invalid section length");
        }

        len -= 2;

        Section section{.marker = marker, .length = len};
        std::vector<uint8_t> data;

        while (len--) {
            data.push_back(reader.Expect1b());
        }

        if (marker == MarkerType::SOS) {
            while (true) {
                uint8_t byte = reader.Expect1b();

                if (byte == kMarkerPrefix) {
                    byte = reader.Expect1b();

                    if (!byte) {
                        data.push_back(kMarkerPrefix);
                    } else {
                        if (SuffixToMarker(byte) != MarkerType::EOI) {
                            throw std::invalid_argument(
                                "Parse: Met unexpected marker in SOS section");
                        }

                        break;
                    }
                } else {
                    data.push_back(byte);
                }
            }

            done = true;
        }

        if (!IsAppnMarker(marker_val)) {
            section.data = std::move(data);
            sections.emplace_back(std::move(section));
        }
    }

    JPEGData rval{.encoding_type = encoding_type, .sections = sections};

#ifdef ENABLE_LOG
    DLOG(INFO) << "Parse: Successfully divided data into sections";
    DLOG(INFO) << ParseResultLogString(rval);
#endif

    return rval;
}

std::string ParseResultLogString(const JPEGData& data) {
    std::stringstream ss;
    ss << std::setfill('0');

    ss << "\nEncoding type: ";

    if (data.encoding_type == EncodingType::BASELINE) {
        ss << "BASELINE";
    } else if (data.encoding_type == EncodingType::PROGRESSIVE) {
        ss << "Progressive";
    } else {
        ss << "Unknown";
    }

    ss << "\n"
       << "Sections number: " << data.sections.size() << "\n";

    for (const auto& section : data.sections) {
        ss << "Marker: " << std::setw(4) << std::hex << MarkerToVal(section.marker) << std::dec
           << "\n";
        ss << "Length: " << section.length << "\n";

        for (size_t i = 0; i < section.data.size(); ++i) {
            if (i != 0 && i % 16 == 0) {
                ss << "\n";
            }

            ss << std::setw(2) << std::hex << static_cast<uint16_t>(section.data[i]) << std::dec
               << " ";
        }

        ss << "\n\n";
    }

    return ss.rdbuf()->str();
}