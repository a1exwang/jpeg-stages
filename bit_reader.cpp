#include "bit_reader.hpp"

int BitReader::read_bit() {
  if (stream_remaining_bits == 0) {
    // std::string has one more byte for '\0' so we are safe
//    octet = data[current_offset++];
    octet = is.get();
    if (octet == 0xFF) {
//      auto b = data[current_offset++];
      auto b = is.get();
      CHECK(b == 0);
    }
//    if (current_offset >= data.size()) {
    if (!is) {
      return -1;
    }
    stream_remaining_bits = 8;
  }
  // Take the MSB(Most significant bit)
  auto ret = (octet & (1 << (stream_remaining_bits-1))) ? 1 : 0;
  stream_remaining_bits--;
  printf("readbit: %d\n", ret);
  return ret;
}

int BitReader::read_nbits_fast(int n, int &output) {
  output = 0;
  int rest = n;
  while (rest > stream_remaining_bits) {
    rest -= stream_remaining_bits;
    output |= (octet >> (8-stream_remaining_bits)) << rest;
//    octet = data[current_offset++];
    octet = is.get();
    if (octet == 0xFF) {
//      auto b = data[current_offset++];
      auto b = is.get();
      CHECK(b == 0);
    }

    // old_octet[0..stream_remaining_bits-1], octet[8..8-n1+1]
//    if (current_offset >= data.size()) {
    if (!is) {
      output >>= rest;
      return n - rest;
    }
    stream_remaining_bits = 8;
  }
  // ASSERT rest <= stream_remaining_bits <= 8
  // Take the MSB(Most significant bit)
  auto mask = (1 << rest) - 1;
  output |= (octet >> (8 - rest)) & mask;
  stream_remaining_bits -= rest;
  return n;
}

int BitReader::read_nbits_safe(int n) {
  int ret = 0;
  for (int i = 0; i < n; i++) {
    ret |= (read_bit() << (n-i-1));
  }
  return ret;
}
