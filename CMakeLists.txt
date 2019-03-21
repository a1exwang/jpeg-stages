cmake_minimum_required(VERSION 3.13)
project(jpeg-stages)

set(CMAKE_CXX_STANDARD 17)
# Currently Arch Linux OpenCV 4.0 CMake is not working, so we manually add a include path
include_directories(/usr/include/opencv4)
include_directories(${CMAKE_SOURCE_DIR})

add_library(jpegdecoder decode.cpp huffman.cpp bit_reader.cpp)

add_executable(jpeg_stages main.cpp)
target_link_libraries(jpeg_stages glog fftw3 opencv_highgui opencv_imgcodecs opencv_core pthread jpegdecoder)

find_package(GTest)
add_executable(tests test/unit/test_bit_reader.cpp)
target_link_libraries(tests glog fftw3 opencv_highgui opencv_imgcodecs opencv_core pthread jpegdecoder gtest ${GTEST_MAIN_LIBRARY})

add_executable(fftw_multi fftw_multi.cpp)
target_link_libraries(fftw_multi fftw3 pthread)
