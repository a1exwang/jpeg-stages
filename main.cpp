#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <thread>
#include <atomic>
#include <glog/logging.h>

#include "decode.hpp"

using namespace std;


int main(int argc, char **argv) {
  if (argc < 3) {
    cerr << "Wrong arguments" << endl;
    return 1;
  }
  FLAGS_logtostderr = true;
  google::InitGoogleLogging(argv[0]);

  string file_path = argv[1];
  string operation = argv[2];
  ifstream ifs(file_path, ios::binary);
  if (!ifs) {
    cerr << "Bad file" << endl;
    return 1;
  }
  string jpeg_data((istreambuf_iterator<char>(ifs)), istreambuf_iterator<char>());


  if (operation == "tp") {
    int test_times = atoi(argv[3]);
    int64_t bs = atoi(argv[4]);
    atomic<int64_t> counter{0};
    atomic<bool> exit{false};
    int64_t last_n = 0;

    thread t([&]() {
      while (true) {
        if (exit) {
          break;
        }
        this_thread::sleep_for(chrono::seconds(1));
        int64_t n = counter;
        cout << (n - last_n) << "files/s" << endl;

        last_n = n;
      }
    });
    for (int i = 0; i < test_times; i++) {
      jst_decode(jpeg_data, bs);
      counter++;
    }
    exit = true;
    t.join();
  } else if (operation == "show") {
    int nbits = atoi(argv[3]);
    cv::Mat mat = jst_decode(jpeg_data, nbits);
    cv::Mat original = cv::imread(file_path);
    cv::Mat roi = cv::Mat(mat.rows, mat.cols*2 + 100, CV_8UC3);
    auto roi1 = roi(cv::Rect(0, 0, original.cols, original.rows));
    original.copyTo(roi1);
    auto roi2 = roi(cv::Rect(mat.cols+100, 0, mat.cols, mat.rows));
    mat.copyTo(roi2);
    cv::imshow("test", roi);
    cv::waitKey(0);
  } else if (operation == "compare") {
    int times = atoi(argv[3]);
    int bs = atoi(argv[4]);
    auto t0 = std::chrono::high_resolution_clock::now();
    cv::Mat original;
    for (int i = 0; i < times; i++) {
      original = cv::imread(file_path);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    int crop_w = 3, crop_h = 3;
    auto t2 = std::chrono::high_resolution_clock::now();
    cv::Mat mat;
    for (int i = 0; i < times; i++) {
      mat = jst_decode(jpeg_data, bs);
    }
    auto t3 = std::chrono::high_resolution_clock::now();
    cout << "showing diff of first " << crop_w << "x" << crop_h << " values" << endl;
    for (int i = 0; i < crop_w; i++) {
      for (int j = 0; j < crop_h; j++) {
        for (int k = 0; k < 3; k++) {
          auto o = original.at<cv::Vec3b>(i, j)[k];
          auto my = mat.at<cv::Vec3b>(i, j)[k];
          cout << (size_t)my << "_" << (size_t)o << " ";
        }
      }
      cout << endl;
    }
    double snr = cv::norm(mat, original, cv::NORM_L2) / (mat.cols*mat.rows*3*255);
    auto d_cv = chrono::duration<double>(t1 - t0).count();
    auto d_my = chrono::duration<double>(t3 - t2).count();
    cout << "time: OpenCV: " << d_cv << "s, my: " << d_my << "s, mine is " << (d_my > d_cv ? d_my/d_cv : d_cv/d_my) << "x " << (d_my>d_cv ? "slower":"faster") << endl;
    double db = 10*log(snr)/log(10);
    cout << "SNR: " << snr << " " << db << "dB" << endl;
    float db_limit = -50;
    if (db > db_limit) {
      cout << "Warning!!!!!!!!!!!!!!!! SNR is greater than " << db_limit << "dB" << endl;
    }
  } else if (operation == "decode") {
    jst_decode(jpeg_data, 2);
  }


  return 0;
}
