#pragma once
#include "bit_reader.hpp"
#include <vector>
#include <map>
#include <cstdint>
#include <optional>
#include <iostream>
#include <glog/logging.h>
#include <queue>

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
    Node(NodeValue value, int64_t index) :value(value), index(index) { }
    NodeValue value;
    std::vector<std::vector<int64_t>> subtrees;
    int64_t children[2]{-1, -1};
    int64_t index;
  };

  explicit HuffmanDecoder(
      std::string str,
      std::map<int64_t, HuffmanTable> dc_dhts,
      std::map<int64_t, HuffmanTable> ac_dhts,
      int64_t bit_batch_size = 4
  );

  std::vector<int> HuffmanDecode(int64_t dc_index, int64_t ac_index, int64_t count_to_read);
  std::vector<int> HuffmanDecodeDFA(int64_t dc_index, int64_t ac_index, int64_t count_to_read);
  uint8_t read_tree(std::vector<Node> &tree) {
#ifndef NDEBUG
    printf("readtree start\n");
#endif
    auto ret = read_tree_batched(tree);
//    auto ret = read_tree_safe(tree);
    read_counter++;
#ifndef NDEBUG
    printf("readtree: %d n=%ld\n", ret, read_counter);
#endif
    return ret;
  }
  int64_t read_counter = 0;
  uint8_t read_tree_safe(std::vector<Node> &tree);
  uint8_t read_tree_batched(std::vector<Node> &tree);
  int64_t get_offset() { return bit_reader.get_offset(); }
private:

  static Node &left_node(std::vector<Node>& tree, const Node &node) {
//    CHECK(node.left >= 0);
    return tree[node.children[0]];
  }
  static Node &right_node(std::vector<Node>& tree, const Node &node) {
//    CHECK(node.right >= 0);
    return tree[node.children[1]];
  }

  static void convert_table_to_tree(
      std::map<int64_t, HuffmanTable> dht,
      std::vector<std::vector<Node>> &trees);

  static void build_subtree(std::vector<std::vector<Node>> &trees, int max_depth);
  static int64_t take_precedence(std::vector<Node> &tree, Node &node, int64_t child_id, int64_t nbits);

  std::vector<std::vector<Node>> dc_trees;
  std::vector<std::vector<Node>> ac_trees;

  BitReader bit_reader;
  std::queue<uint8_t> cached_values;
  int64_t bit_batch_size;
};

