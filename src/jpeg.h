#pragma once
#include <vector>
#include <string>
#include <cstdint>

enum class EncodingType { BASELINE, PROGRESSIVE };

enum class MarkerType : uint8_t {
    SOI = 0xD8,
    SOF0 = 0xC0,
    SOF2 = 0xC2,
    DHT = 0xC4,
    DQT = 0xDB,
    SOS = 0xDA,
    COM = 0xFE,
    EOI = 0xD9
};

constexpr uint8_t kMarkers[] = {
    static_cast<uint8_t>(MarkerType::SOI),  static_cast<uint8_t>(MarkerType::SOF0),
    static_cast<uint8_t>(MarkerType::SOF2), static_cast<uint8_t>(MarkerType::DHT),
    static_cast<uint8_t>(MarkerType::DQT),  static_cast<uint8_t>(MarkerType::SOS),
    static_cast<uint8_t>(MarkerType::COM),  static_cast<uint8_t>(MarkerType::EOI)};

const uint8_t kMarkerPrefix = 0xFF;

uint16_t MarkerToVal(MarkerType type);
MarkerType ValToMarker(uint16_t val);
uint8_t MarkerToSuffix(MarkerType type);
MarkerType SuffixToMarker(uint8_t suf);
bool IsMarker(uint16_t val);
bool IsAppnMarker(uint16_t val);

struct Section {
    MarkerType marker;
    size_t length;
    std::vector<uint8_t> data;
};

struct JPEGData {
    EncodingType encoding_type;
    std::vector<Section> sections;
};
