#include "stream_reader.h"
#include <cstdint>
#include <stdexcept>

StreamReader::StreamReader(std::istream& in) : in_(in) {
}

bool StreamReader::MustEqual2b(uint16_t x) {
    char buf[2];
    in_.read(buf, 2);

    if (!in_) {
        return false;
    }

    for (size_t i = 0; i < 2; ++i) {
        if (*reinterpret_cast<uint8_t*>(buf + i) != ((x >> ((1 - i) * 8)) & 0xFF)) {
            return false;
        }
    }

    return true;
}

uint16_t StreamReader::Expect2b() {
    char buf[2];
    in_.read(buf, 2);

    if (!in_) {
        throw std::runtime_error("StreamReader::Expect2b: Couldn't read enough bytes");
    }

    return (*reinterpret_cast<uint8_t*>(buf) << 8) | *reinterpret_cast<uint8_t*>(buf + 1);
}

bool StreamReader::Read1b(uint8_t* byte) {
    char c = in_.get();

    if (!in_) {
        return false;
    }

    *byte = *reinterpret_cast<uint8_t*>(&c);

    return true;
}

uint8_t StreamReader::Expect1b() {
    char c = in_.get();

    if (!in_) {
        throw std::runtime_error("StreamReader::Expect1b: No byte while one was expected");
    }

    return *reinterpret_cast<uint8_t*>(&c);
}