#include "bit_reader.hpp"

int BitReader::read_nbits_safe(int n, int &output) {
  output = 0;
  for (int i = 0; i < n; i++) {
    output |= (read_bit() << (n-i-1));
  }
  return n;
}
