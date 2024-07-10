#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <filesystem>

#include <decoder.h>
#include <png_encoder.hpp>

Image DecodeJPEGFile(std::string path) {
    std::ifstream in(path);
    if (!in.good()) {
        throw std::runtime_error(std::string("failed to open file: ") + path);
    }
    return Decode(in);
}

int main(int argc, char** args) {
    if (argc != 2) {
        std::cerr << "Usage: ./jpeg-decoder <FILE>\n";
        return 0;
    }

    try {
        std::string in(args[1]);
        std::string out(std::filesystem::path(in).replace_extension("png"));

        Image image = DecodeJPEGFile(in);
        WritePng(out, image);
    } catch (std::exception& e) {
        std::cout << "failed to decode: " << e.what() << "\n";
    } catch (...) {
        std::cout << "fatal error\n";
        return 1;
    }
}
