#pragma once
#include <vector>
#include <map>
#include <cstdint>
#include <iostream>

class HuffmanTable {
public:
  int64_t id;
  bool is_ac;
  std::vector<size_t> counts;
  std::vector<std::vector<int8_t>> table;
  std::vector<uint8_t> codes;
};

class HuffmanDecoder {
public:
  explicit HuffmanDecoder(
      std::istream &is,
      std::map<int64_t, HuffmanTable> dc_dhts,
      std::map<int64_t, HuffmanTable> ac_dhts
  );

  std::vector<int> HuffmanDecode(int64_t dc_index, int64_t ac_index, int64_t count_to_read);
  int read_bit();

private:
  int octet = 0;
  int remaining_bits = 0;
  std::istream &is;
  std::vector<std::vector<int64_t>> dc_trees, ac_trees;
};

