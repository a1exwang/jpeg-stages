cmake_minimum_required(VERSION 3.13)
project(jpeg-stages)

set(CMAKE_CXX_STANDARD 17)
# Currently Arch Linux OpenCV 4.0 CMake is not working, so we manually add a include path
include_directories(/usr/include/opencv4)

add_library(jpegdecoder decode.cpp huffman.cpp)

add_executable(jpeg_stages main.cpp)
target_link_libraries(jpeg_stages glog fftw3 opencv_highgui opencv_imgcodecs opencv_core pthread jpegdecoder)

add_executable(fftw_multi fftw_multi.cpp)
target_link_libraries(fftw_multi fftw3 pthread)
