#include "bit_reader.hpp"

int BitReader::read_bit() {
  if (stream_remaining_bits == 0) {
    // std::string has one more byte for '\0' so we are safe
    octet = read_byte();
    if (octet == 0xFF) {
      auto b = read_byte();
      CHECK(b == 0);
    }
    if (next_offset > data.size()) {
      return -1;
    }
    stream_remaining_bits = 8;
  }
  // Take the MSB(Most significant bit)
  auto ret = (octet & (1 << (stream_remaining_bits-1))) ? 1 : 0;
  stream_remaining_bits--;
#ifndef NDEBUG
  printf("offset: %ld, bitoff: %d, readbit: %d\n", next_offset-1, 7-stream_remaining_bits, ret);
#endif
  return ret;
}

int BitReader::read_nbits_fast(int n, int &output) {
  output = 0;
  int rest = n;
  while (rest >= stream_remaining_bits) {
    if (stream_remaining_bits == 0) {
      octet = read_byte();
      if (octet == 0xFF) {
        auto b = read_byte();
        CHECK(b == 0);
      }
    }

    rest -= stream_remaining_bits;
    auto mask = (1 << stream_remaining_bits) - 1;
    auto tmp = ((octet >> (8-stream_remaining_bits)) & mask);
    if (stream_remaining_bits > 0) {
      printf("read_nbits: %x, offset=%ld, bitoff=%d, n=%d\n", tmp, next_offset-1, 7-stream_remaining_bits, stream_remaining_bits);
    }
    output |= (tmp << rest);
    octet = read_byte();
    if (octet == 0xFF) {
      auto b = read_byte();
      CHECK(b == 0);
    }

    // old_octet[0..stream_remaining_bits-1], octet[8..8-n1+1]
    if (next_offset > data.size()) {
      output >>= rest;
      return n - rest;
    }
    stream_remaining_bits = 8;
  }
  // ASSERT rest <= stream_remaining_bits <= 8
  // Take the MSB(Most significant bit)
  auto mask = (1 << rest) - 1;
  auto tmp = (octet >> (stream_remaining_bits - rest)) & mask;
  if (rest > 0) {
    printf("read_nbits: %x, offset=%ld, bitoff=%d, n=%d\n", tmp, next_offset-1, 8-stream_remaining_bits, rest);
  }
  output |= tmp;
  stream_remaining_bits -= rest;
  return n;
}

int BitReader::read_nbits_safe(int n, int &output) {
  output = 0;
  for (int i = 0; i < n; i++) {
    output |= (read_bit() << (n-i-1));
  }
  return n;
}
