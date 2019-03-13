#pragma once

#include <opencv2/opencv.hpp>
#include <string>
cv::Mat jst_decode(const std::string &jpeg_data, bool verbose = false);
