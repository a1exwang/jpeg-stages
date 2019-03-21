#pragma once
#include <iostream>
#include <glog/logging.h>

class BitReader {
public:
  explicit BitReader(std::istream &is) :is(is) { }
  int read_bit();
  int read_nbits(int n) { return read_nbits_fast(n); }
  int read_nbits_fast(int n);
  int read_nbits_safe(int n);
private:
  int octet = 0;
  int stream_remaining_bits = 0;
  std::istream &is;
};