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

enum class NodeValue :int64_t {
  NonExisting = -1,
  NonExistingSubtree = -2,
};

class HuffmanDecoder {
public:
  struct Node {
    Node(NodeValue value) :value(value) { }
    NodeValue value;
    NodeValue subtree_l1[2];
    NodeValue subtree_l2[4];
    int64_t left{-1};
    int64_t right{-1};
  };

  explicit HuffmanDecoder(
      std::string str,
      std::map<int64_t, HuffmanTable> dc_dhts,
      std::map<int64_t, HuffmanTable> ac_dhts
  );

  std::vector<int> HuffmanDecode(int64_t dc_index, int64_t ac_index, int64_t count_to_read);
  std::vector<int> HuffmanDecodeDFA(int64_t dc_index, int64_t ac_index, int64_t count_to_read);
  uint8_t read_tree(std::vector<Node> &tree);
  uint8_t read_tree_batched(std::vector<Node> &tree);
  int64_t get_offset() { return bit_reader.get_offset(); }
private:

  static Node &left_node(std::vector<Node>& tree, const Node &node) {
    CHECK(node.left >= 0);
    return tree[node.left];
  }
  static Node &right_node(std::vector<Node>& tree, const Node &node) {
    CHECK(node.right >= 0);
    return tree[node.right];
  }

  static void convert_table_to_tree(
      std::map<int64_t, HuffmanTable> dht,
      std::vector<std::vector<Node>> &trees);

  static void build_subtree(std::vector<std::vector<Node>> &trees);

  std::vector<std::vector<Node>> dc_trees;
  std::vector<std::vector<Node>> ac_trees;

  BitReader bit_reader;
};

