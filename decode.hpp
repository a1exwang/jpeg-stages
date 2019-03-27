#pragma once

#include <opencv2/opencv.hpp>
#include <string>
cv::Mat jst_decode(const std::string &jpeg_data, int64_t decode_batch_size, bool verbose = false);
