cmake_minimum_required(VERSION 3.13)

project(JpegDecoder)

option(BUILD_TESTS OFF)

set(CMAKE_MODULE_PATH              "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_EXPORT_COMPILE_COMMANDS  ON)

include(cmake/BuildFlags.cmake)

if(BUILD_TESTS)
  include(cmake/TestSolution.cmake)
endif()

include(FetchContent)
include(cmake/FindGlog.cmake)

if(BUILD_TESTS)
  find_package(Catch REQUIRED)
endif()

include_directories(util)

add_library(decoder_baseline
        src/huffman.cpp
        src/fft.cpp
        src/parser.cpp
        src/stream_reader.cpp
        src/jpeg.cpp
        src/decoder.cpp
)

target_include_directories(decoder_baseline PRIVATE src)

find_package(PNG)
find_package(FFTW)

if(BUILD_TESTS)
  find_package(JPEG)
endif()

message(STATUS "Path to PNG includes: ${PNG_INCLUDE_DIRS}")
message(STATUS "Path to PNG library: ${PNG_LIBRARY}")
message(STATUS "Path to JPEG includes: ${JPEG_INCLUDE_DIRS}")
message(STATUS "Path to JPEG library: ${JPEG_LIBRARIES}")
message(STATUS "Path to FFTW includes: ${FFTW_INCLUDES}")
message(STATUS "Path to FFTW library: ${FFTW_LIBRARIES}")

set(DECODER_UTIL_FILES
    decoder-utils/logger_init.cpp
    decoder-utils/libjpg_reader.cpp
    decoder-utils/png_encoder.cpp 
    decoder-utils/test_commons.cpp
)


if(BUILD_TESTS)
    add_catch(test_decoder_baseline
        src/tests/test_huffman.cpp
        src/tests/test_fft.cpp
        src/tests/test_baseline.cpp
        ${DECODER_UTIL_FILES}
    )

    target_compile_definitions(test_decoder_baseline PUBLIC HSE_TASK_DIR="${CMAKE_CURRENT_SOURCE_DIR}/")
endif()

set(DECODER_UTILS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/decoder-utils)
set(DECODER_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

function(link_decoder_deps TARGET)
    target_include_directories(${TARGET} PUBLIC
            ${PNG_INCLUDE_DIRS}
            ${FFTW_INCLUDES}
            ${DECODER_UTILS_DIR}
            ${DECODER_INCLUDE_DIR})

    target_link_libraries(${TARGET} PUBLIC
            ${FFTW_LIBRARIES}
            ${PNG_LIBRARY}
            glog::glog)
    if(BUILD_TESTS)
        target_include_directories(${TARGET} PUBLIC ${JPEG_INCLUDE_DIRS})
        target_link_libraries(${TARGET} PUBLIC ${JPEG_LIBRARIES})
    endif()

    get_target_property(GLOG_INCLUDES glog::glog INCLUDE_DIRECTORIES)
    target_include_directories(${TARGET} SYSTEM PUBLIC ${GLOG_INCLUDES})

endfunction()

link_decoder_deps(decoder_baseline)

if(BUILD_TESTS)
    target_link_libraries(test_decoder_baseline decoder_baseline)
    target_include_directories(test_decoder_baseline PUBLIC ${PROJECT_SOURCE_DIR}/src)
endif()

add_executable(jpeg-decoder src/main.cpp decoder-utils/png_encoder.cpp)
target_link_libraries(jpeg-decoder PUBLIC decoder_baseline)
target_include_directories(jpeg-decoder PUBLIC ${PROJECT_SOURCE_DIR}/src)
