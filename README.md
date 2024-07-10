# JPEG Decoder

Just as the title suggests, this small piece of software can decode `.jpeg` or `.jpg` files for you.

Once built, you can use the binary like this:

    ./jpeg-decoder <FILE>

It will output the decoded image as a `.png` file.

## Building from Source

To build the main binary, you need a few dependencies: libpng, glog and fftw3. 
Additionaly, if you wish to build tests you need to specify 
`-DBUILD_TESTS=ON` CMake option and have libjpeg as a dependency.

Here's an example of how you can build the main executable:

    cmake -DCMAKE_BUILD_TYPE=Release -S . -B build
    cd build && cmake --build .
