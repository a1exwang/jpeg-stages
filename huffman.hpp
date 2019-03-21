#pragma once
#include "bit_reader.hpp"
#include <vector>
#include <map>
#include <cstdint>
#include <iostream>
#include <glog/logging.h>

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
  struct Node {
    Node(int64_t value) :value(value) { }
    int64_t value;
  };

  explicit HuffmanDecoder(
      std::istream &is,
      std::map<int64_t, HuffmanTable> dc_dhts,
      std::map<int64_t, HuffmanTable> ac_dhts
  );

  std::vector<int> HuffmanDecode(int64_t dc_index, int64_t ac_index, int64_t count_to_read);
  uint8_t read_tree(std::vector<Node> &tree);

private:
  static void convert_table_to_tree(
      std::map<int64_t, HuffmanTable> dht,
      std::vector<std::vector<Node>> &trees);

  std::vector<std::vector<Node>> dc_trees;
  std::vector<std::vector<Node>> ac_trees;

  BitReader bit_reader;
};

