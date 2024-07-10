#include "jpeg.h"
#include <array>

namespace {
constexpr std::array<bool, 256> MarkMarkers() {
    std::array<bool, 256> marked = {};

    for (size_t i = 0; i < sizeof(kMarkers); ++i) {
        marked[static_cast<uint8_t>(kMarkers[i])] = true;
    }

    return marked;
}

constexpr std::array<bool, 256> kMarked = MarkMarkers();
}  // namespace

uint16_t MarkerToVal(MarkerType type) {
    return (static_cast<uint16_t>(kMarkerPrefix) << 8) | static_cast<uint8_t>(type);
}

MarkerType ValToMarker(uint16_t val) {
    return static_cast<MarkerType>(val & 0xff);
}

uint8_t MarkerToSuffix(MarkerType type) {
    return static_cast<uint8_t>(type);
}

MarkerType SuffixToMarker(uint8_t suf) {
    return static_cast<MarkerType>(suf);
}

bool IsMarker(uint16_t val) {
    return IsAppnMarker(val) || (kMarked[val & 0xff] && ((val >> 8) == 0xff));
}

bool IsAppnMarker(uint16_t val) {
    return ((val ^ 0xffe0) & 0xfff0) == 0;
}
