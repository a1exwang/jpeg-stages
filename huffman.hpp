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
    Node() :value(NodeValue::NonExisting), index(-1) { }
    Node(NodeValue value, int64_t index) :value(value), index(index) { }
    NodeValue value;
    std::vector<std::vector<Node*>> subtrees;
    Node *children[2]{nullptr, nullptr};
    int64_t index;

    Node(const NodeValue &) = delete;
    Node &operator=(const NodeValue &) = delete;
  };

  explicit HuffmanDecoder(
      std::string str,
      std::map<int64_t, HuffmanTable> dc_dhts,
      std::map<int64_t, HuffmanTable> ac_dhts,
      int64_t bit_batch_size = 4
  );

  std::vector<int> HuffmanDecode(int64_t dc_index, int64_t ac_index, int64_t count_to_read);
  uint8_t read_tree(std::vector<Node> &tree, int &zzk, bool is_dc) {
#ifndef NDEBUG
    printf("readtree {\n");
#endif
    auto ret = read_tree_batched(tree, zzk, is_dc);
//    auto ret = read_tree_safe(tree);
    read_counter++;
#ifndef NDEBUG
    printf("} %d n=%ld\n", ret, read_counter);
#endif
    return ret;
  }
  int64_t read_counter = 0;
  uint8_t read_tree_safe(std::vector<Node> &tree);
  uint8_t read_tree_batched(std::vector<Node> &tree, int &zzk, bool is_dc);
  bool read_tree_nbits(Node **node, int &batch_size, int bits, uint8_t &value);
  uint8_t read_tree_fallback(std::vector<Node> &tree, Node *current, int nread, int bits) {
    // fallback

    Node *node = current;
    for (int i = 0; i < nread; i++) {
      int bit = (bits >> (nread - i - 1)) & 1;
      // 0 for left, 1 for right
      node = node->children[bit];
      NodeValue value = node->value;
      if (value != NodeValue::NonExisting) {
        // return remaining bits
        bit_reader.return_nbits(nread - i - 1);
#ifndef NDEBUG
        printf("readtree: return bits: %d\n", nread-i-1);
#endif
        return (uint8_t)value;
      }
    }
    LOG(FATAL) << "Unreachable code";
    return 0;
  }
  int64_t get_offset() { return bit_reader.get_offset(); }
private:
  static void convert_table_to_tree(
      std::map<int64_t, HuffmanTable> dht,
      std::vector<std::vector<Node>> &trees);

  static void build_subtree(std::vector<std::vector<Node>> &trees, int max_depth);
  static Node *take_precedence(std::vector<Node> &tree, Node &node, int64_t child_id, int64_t nbits);

  std::vector<std::vector<Node>> dc_trees;
  std::vector<std::vector<Node>> ac_trees;

  BitReader bit_reader;
  std::queue<uint8_t> cached_values;
  int64_t bit_batch_size;
};

