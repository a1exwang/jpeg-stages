#pragma once
#include <iostream>
#include <glog/logging.h>

class BitReader {
public:
  explicit BitReader(std::istream &is) :is(is) { }
//  explicit BitReader(std::string str) :data(std::move(str)) { }
  int read_bit();
  int read_nbits(int n, int &output) { return read_nbits_fast(n, output); }
  int read_nbits_fast(int n, int &output);
  int read_nbits_safe(int n);
  int64_t get_offset() { return current_offset; }
private:
  int octet = 0;
  int stream_remaining_bits = 0;
  std::istream &is;
//  std::string data;
  int64_t current_offset{0};
};