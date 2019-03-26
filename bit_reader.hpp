#pragma once
#include <iostream>
#include <glog/logging.h>
#include <cmath>

class BitReader {
public:
  explicit BitReader(std::string str) :data(std::move(str)) { }
  int read_bit();
  int read_nbits(int n, int &output) { return read_nbits_fast(n, output); }
  int read_nbits_fast(int n, int &output);
  int read_nbits_safe(int n, int &output);
  void return_bit() {
    if (stream_remaining_bits < 7) {
      stream_remaining_bits++;
      return;
    }
    CHECK(stream_remaining_bits == 7);

    next_offset--;
    if (octet == 0xff) {
      next_offset--;
    }
    octet = (int)(uint32_t)(uint8_t)data[next_offset-1];
    stream_remaining_bits = 0;
  }
  void return_nbits(int n);
  int64_t get_offset() { return next_offset; }
private:
  int read_byte() {
    auto ret = (int)(uint32_t)(uint8_t)data[next_offset++];
#ifndef NDEBUG
    printf("read_byte: %x\n", ret);
#endif
    return ret;
  }
  int octet = 0;
  int stream_remaining_bits = 0;
  std::string data;
  int64_t next_offset{0};

};