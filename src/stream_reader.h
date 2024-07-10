#pragma once
#include <istream>
#include <array>

class StreamReader {
public:
    explicit StreamReader(std::istream& in);
    StreamReader(StreamReader& o) = delete;
    StreamReader(StreamReader&& o) = delete;
    StreamReader& operator=(StreamReader& o) = delete;
    StreamReader& operator=(StreamReader&& o) = delete;
    ~StreamReader() = default;

    bool MustEqual2b(uint16_t x);
    uint16_t Expect2b();
    bool Read1b(uint8_t* byte);
    uint8_t Expect1b();

private:
    std::istream& in_;
};