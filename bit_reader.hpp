#pragma once
#include <iostream>
#include <glog/logging.h>
#include <cmath>

class BitReader {
public:
  explicit BitReader(std::string str) :data(std::move(str)) { }
  int read_bit() {
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
  int read_nbits(int n, int &output) { return read_nbits_fast(n, output); }
  static void print_bits(std::ostream &os, int bits, int n) {
    CHECK(n > 0);
    int mask = 1 << (n-1);
    for (int i = 0; i < n; i ++) {
      os << (((bits & mask) == 0) ? "0" : "1");
      mask >>= 1;
    }
  }
  int read_nbits_fast(int n, int &output) {
#ifndef NDEBUG
    int64_t original_offset = next_offset-1;
    int original_bit_offset = 8-stream_remaining_bits;

    std::cout << "  read_nbits:"
              << " offset0=" << original_offset
              << " bitoff0=" << original_bit_offset << " {" << std::endl;
#endif

    output = 0;
    int rest = n;
    while (rest >= stream_remaining_bits) {
      rest -= stream_remaining_bits;
      output |= (octet & ((1 << stream_remaining_bits) - 1)) << rest;
      if (rest == 0) {
        stream_remaining_bits = 0;
        return n;
      }

      octet = read_byte();
      if (octet == 0xFF) {
        auto b = read_byte();
        CHECK(b == 0);
      }

      // old_octet[0..stream_remaining_bits-1], octet[8..8-n1+1]
      if (next_offset > data.size()) {
        output >>= rest;
#ifndef NDEBUG
        if (n - rest > 0) {
          std::stringstream ss;
          print_bits(ss, output, n-rest);
          std::cout << "  } n=" << n-rest
                    << " value=0b" << ss.str()
                    << "(0x" << std::hex << output << std::dec << ")" << std::endl;
        }
#endif
        return n - rest;
      }
      stream_remaining_bits = 8;
    }
    // ASSERT rest <= stream_remaining_bits <= 8
    // Take the MSB(Most significant bit)
    auto tmp = (octet >> (stream_remaining_bits - rest)) & ((1 << rest) - 1);
    output |= tmp;
    stream_remaining_bits -= rest;

#ifndef NDEBUG
    if (n > 0) {
      std::stringstream ss;
      print_bits(ss, output, n);
      std::cout << "  }"
                << " n=" << n
                << " value=0b" << ss.str()
                << "(0x" << std::hex << output << std::dec << ")" << std::endl;
    }
#endif
    return n;
  }
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
  void return_nbits(int n) {
#ifndef NDEBUG
    std::cout << "  return_nbits: " << n << std::endl;
#endif
    if (n > (8-stream_remaining_bits)) {
      n -= (8-stream_remaining_bits);
      next_offset--;
      stream_remaining_bits = 0;
      if (next_offset >= 1 && data[next_offset - 1] == -1) {
        next_offset--;
      }
    }

    stream_remaining_bits += n % 8;
    next_offset -= int64_t(n / 8);
    if (int64_t(n / 8) != 0) {
      if (next_offset >= 1 && data.at(next_offset - 1) == -1) {
        next_offset--;
      }
    }
    if (next_offset >= 2 && data.at(next_offset - 2) == -1) {
      octet = data[next_offset-2];
    } else {
      octet = data[next_offset-1];
    }
//    CHECK(next_offset >= 0);
  }
  int64_t get_offset() { return next_offset; }
private:
  int read_byte() {
    auto ret = (int)(uint32_t)(uint8_t)data[next_offset++];
#ifndef NDEBUG
    printf("    read_byte: %x\n", ret);
#endif
    return ret;
  }
  int octet = 0;
  int stream_remaining_bits = 0;
  std::string data;
  int64_t next_offset{0};
};