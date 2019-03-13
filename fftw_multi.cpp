#include <fftw3.h>
#include <iostream>
#include <atomic>
#include <string>
#include <vector>
#include <thread>

constexpr int64_t width = 8;
using namespace std;

int main(int argc, char **argv) {
  if (argc < 2) {
    cerr << "invalid argument" << endl;
    return 1;
  }

  int64_t batch_size = atoi(argv[1]);
  int64_t width = atoi(argv[2]);
  vector<double> input(width*width*batch_size, 0), output(width*width*batch_size, 0);
  vector<int> ns(2, width);
  vector<fftw_r2r_kind> kinds(batch_size, FFTW_REDFT01);

  atomic<int64_t> total_size{0};
  int64_t last_n = 0;
  thread([&]() {
    int quit = 0;
    while (quit < 5) {
      this_thread::sleep_for(chrono::seconds(1));
      int64_t n = total_size;
      cout << "tp: " << (n - last_n)/1e6 << "MB/s" << endl;

      last_n = n;
      quit++;
    }
    exit(0);
  }).detach();


  while (true) {

    auto plan = fftw_plan_many_r2r(
        2, ns.data(), batch_size, input.data(),
        nullptr, 0, 0, output.data(), nullptr, 0, 0, kinds.data(), 0);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    total_size += input.size();
  }
}