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
      jst_decode(jpeg_data);
      counter++;
    }
    exit = true;
    t.join();
  } else if (operation == "show") {
    cv::Mat mat = jst_decode(jpeg_data);
    cv::imshow("test", mat);
    cv::waitKey(0);
  } else if (operation == "compare") {
    auto t0 = std::chrono::high_resolution_clock::now();
    auto original = cv::imread(file_path);
    auto t1 = std::chrono::high_resolution_clock::now();
    int crop_w = 10, crop_h = 10;
    auto t2 = std::chrono::high_resolution_clock::now();
    cv::Mat mat = jst_decode(jpeg_data);
    auto t3 = std::chrono::high_resolution_clock::now();
    cout << "showing diff of first " << crop_w << "x" << crop_h << " values" << endl;
    for (int i = 0; i < crop_w; i++) {
      for (int j = 0; j < crop_h; j++) {
        auto o = original.at<cv::Vec3b>(i, j)[0];
        auto my = mat.at<cv::Vec3b>(i, j)[0];
        cout << (size_t)my << "_" << (size_t)o << " ";
      }
      cout << endl;
    }
    auto d_cv = chrono::duration<double>(t1 - t0).count();
    auto d_my = chrono::duration<double>(t3 - t2).count();
    cout << "time: OpenCV: " << d_cv << "s, my: " << d_my << "s, mine is " << (d_my > d_cv ? d_my/d_cv : d_cv/d_my) << "x " << (d_my>d_cv ? "slower":"faster") << endl;
  } else if (operation == "decode") {
    jst_decode(jpeg_data);
  }


  return 0;
}
