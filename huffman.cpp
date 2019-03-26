#include <vector>
#include "huffman.hpp"
#include <string>
#include <iostream>
#include "decode.hpp"
#include <glog/logging.h>

using namespace std;


constexpr int64_t MaximumHuffmanDepth = 16;

static inline int extend(int value, int nBits) {
  // value < 2^nBits
  if (value >> (nBits - 1) == 0) {
    value -= (1 << nBits) - 1;
  }
  return value;
}

std::vector<int> HuffmanDecoder::HuffmanDecode(int64_t dc_index, int64_t ac_index, int64_t count_to_read) {
  vector<int> ret(count_to_read, 0);
  int64_t count_read = 0;
  // Finish building huffman tree
  auto t = read_tree(dc_trees[dc_index]);
  int diff = 0;
  if (t == 0) {
    diff = 0;
  } else {
    bit_reader.read_nbits(t, diff);
    diff = extend(diff, t);
  }

  ret.at(count_read++) = diff;

  while (count_read < count_to_read) {
    auto ac_value = read_tree(ac_trees.at(ac_index));
    int r = (ac_value >> 4) & 0xF;
    int ssss = (ac_value & 0xF);
    if (ssss == 0) {
      if (r == 15) {
        count_read += 16;
        continue;
      } else {
        break;
      }
    }
    // ssss > 0
    count_read += r;

    // decode_zz k
    int zzk = -1;
    bit_reader.read_nbits(ssss, zzk);
    CHECK(count_read < count_to_read);
    ret.at(count_read++) = extend(zzk, ssss);
  }
  return std::move(ret);
}

HuffmanDecoder::HuffmanDecoder(
    std::string str,
    std::map<int64_t, HuffmanTable> dc_dhts,
    std::map<int64_t, HuffmanTable> ac_dhts,
    int64_t bit_batch_size
) :bit_reader(move(str)), dc_trees(dc_dhts.size()), ac_trees(ac_dhts.size()), bit_batch_size(bit_batch_size) {

  convert_table_to_tree(dc_dhts, dc_trees);
  convert_table_to_tree(ac_dhts, ac_trees);
  build_subtree(dc_trees, bit_batch_size);
  build_subtree(ac_trees, bit_batch_size);
}

uint8_t HuffmanDecoder::read_tree_safe(std::vector<Node> &tree) {
  int64_t current = 0;
  while (true) {
    int bit = bit_reader.read_bit();
    if (bit < 0) {
      cerr << "Invalid entropy stream, unexpected EOF" << endl;
      abort();
    }
    // 0 for left, 1 for right
    NodeValue value = tree[tree[current].children[bit]].value;
    current = tree[current].children[bit];
    if (value != NodeValue::NonExisting) {
      // yield data
      return (uint8_t)value;
    }
  }
}

void HuffmanDecoder::convert_table_to_tree(
    std::map<int64_t, HuffmanTable> dht,
    vector<vector<HuffmanDecoder::Node>> &trees
    ) {


  for (int i = 0; i < dht.size(); i++) {
    auto &ac_values = dht[i].codes;
    auto &ac_counts = dht[i].counts;
    auto &tree = trees[i];
    auto it = ac_values.begin();
    size_t skip = 0;
    tree.emplace_back(NodeValue::NonExisting, tree.size());
    size_t last_mid_node_start = 0;
    for (int depth = 0; depth < ac_counts.size(); depth++) {
      skip = (skip + (depth == 0 ? 0 : ac_counts[depth-1]))* 2;
      for (int j = 0; j < /*width*/ ac_counts[depth]; j++) {
        tree.emplace_back(NodeValue(*it++), tree.size());
        // <depth+1, skip+j>

        // if last_mid_node left is empty
        //   put in left
        // else
        //   put in right
        //   last_mid_node_start++;
        if (tree[last_mid_node_start].children[0] == -1) {
          tree[last_mid_node_start].children[0] = tree.size()-1;
        } else {
          tree[last_mid_node_start].children[1] = tree.size()-1;
          last_mid_node_start++;
        }
      }
      auto mid_node_start = tree.size();
      for (size_t j = skip + ac_counts[depth]; j < (1UL << (depth+1)); j++) {
        if (tree[last_mid_node_start].children[0] == -1) {
          tree[last_mid_node_start].children[0] = tree.size();
        } else {
          tree[last_mid_node_start].children[1] = tree.size();
          last_mid_node_start++;
        }
        tree.emplace_back(NodeValue::NonExisting, tree.size());
        // if last_mid_node left is empty
        //   put in left
        // else
        //   put in right
        //   last_mid_node_start++;
      }
      last_mid_node_start = mid_node_start;
    }
  }

}

void HuffmanDecoder::build_subtree(vector<vector<HuffmanDecoder::Node>> &trees, int max_depth) {
  for (auto &tree : trees) {
    for (auto &node : tree) {
      for (int d = 1; d <= max_depth; d++) {
        auto &subtree = node.subtrees.emplace_back(1 << d);
        for (size_t i = 0; i < subtree.size(); i++) {
          subtree[i] = take_precedence(tree, node, i, d);
        }
      }
    }
  }
}

uint8_t HuffmanDecoder::read_tree_batched(vector<HuffmanDecoder::Node> &tree) {
  int64_t current = 0;
  if (!cached_values.empty()) {
    auto ret = cached_values.front();
    cached_values.pop();
    return ret;
  }
  while (true) {
    int bits;
    int nread = bit_reader.read_nbits(bit_batch_size, bits);
//    if (bits < 0) {
//      cerr << "Invalid entropy stream, unexpected EOF" << endl;
//      abort();
//    }
    // 0 for left, 1 for right
//    CHECK(nread == bit_batch_size);
    auto next_node = tree[current].subtrees[bit_batch_size-1].at(bits);
    if (next_node == -1) {
      // fallback
      for (int i = 0; i < nread; i++) {
        int bit = (bits >> (nread - i - 1)) & 1;
        // 0 for left, 1 for right
        auto &node = tree[current];
        current = node.children[bit];
        NodeValue value = tree[current].value;
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

    } else {
      if (int64_t(tree[next_node].value) >= 0) {
        return uint8_t(tree[next_node].value);
      } else {
        current = next_node;
      }
    }
  }
}

int64_t HuffmanDecoder::take_precedence(
    vector<HuffmanDecoder::Node> &tree,
    HuffmanDecoder::Node &node,
    int64_t child_id,
    int64_t nbits) {

  CHECK(nbits <= 63 && nbits > 0);
  Node *current_node = &node;
  for (int64_t i = nbits-1; i >= 0; i--) {
    int64_t current_bit = (child_id >> i) & 0x1;
    auto next_index = current_node->children[current_bit];
//    CHECK(next_index >= 0);
    if (next_index < 0) {
      return -1;
    }

    current_node = &tree[next_index];
  }
  return current_node->index;

}

