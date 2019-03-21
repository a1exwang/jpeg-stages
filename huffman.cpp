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
    ret.at(count_read++) = extend(zzk, ssss);
  }
  return std::move(ret);
}

HuffmanDecoder::HuffmanDecoder(
    std::string str,
    std::map<int64_t, HuffmanTable> dc_dhts,
    std::map<int64_t, HuffmanTable> ac_dhts
) :bit_reader(move(str)), dc_trees(dc_dhts.size()), ac_trees(ac_dhts.size()) {

  convert_table_to_tree(dc_dhts, dc_trees);
  convert_table_to_tree(ac_dhts, ac_trees);
  build_subtree(dc_trees);
  build_subtree(ac_trees);
}

uint8_t HuffmanDecoder::read_tree(std::vector<Node> &tree) {
  int64_t current = 0;
  while (true) {
    int bit = bit_reader.read_bit();
    if (bit < 0) {
      cerr << "Invalid entropy stream, unexpected EOF" << endl;
      abort();
    }
    // 0 for left, 1 for right
    NodeValue value = tree[tree[current].subtree_l1[bit]].value;
    current = tree[current].subtree_l1[bit];
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
        if (tree[last_mid_node_start].subtree_l1[0] == -1) {
          tree[last_mid_node_start].subtree_l1[0] = tree.size()-1;
        } else {
          tree[last_mid_node_start].subtree_l1[1] = tree.size()-1;
          last_mid_node_start++;
        }
      }
      auto mid_node_start = tree.size();
      for (size_t j = skip + ac_counts[depth]; j < (1UL << (depth+1)); j++) {
        if (tree[last_mid_node_start].subtree_l1[0] == -1) {
          tree[last_mid_node_start].subtree_l1[0] = tree.size();
        } else {
          tree[last_mid_node_start].subtree_l1[1] = tree.size();
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

void HuffmanDecoder::build_subtree(vector<vector<HuffmanDecoder::Node>> &trees) {
  for (auto &tree : trees) {
    for (int i = 0; i < tree.size(); i++) {
      if (tree[i].subtree_l1[0] != -1) {
        tree[i].subtree_l2[0b00] = left_node(tree, left_node(tree, tree[i])).index;
        tree[i].subtree_l2[0b10] = right_node(tree, left_node(tree, tree[i])).index;
      }
      if (tree[i].subtree_l1[1] != -1) {
        tree[i].subtree_l2[0b01] = left_node(tree, right_node(tree, tree[i])).index;
        tree[i].subtree_l2[0b11] = right_node(tree, right_node(tree, tree[i])).index;
      }
    }
  }

}

uint8_t HuffmanDecoder::read_tree_batched(vector<HuffmanDecoder::Node> &tree) {
  int64_t current = 0;
  while (true) {
    int bits;
    int nread = bit_reader.read_nbits(2, bits);
    if (bits < 0) {
      cerr << "Invalid entropy stream, unexpected EOF" << endl;
      abort();
    }
    // 0 for left, 1 for right
    NodeValue value;
//    tree[current].subtree_l1
    if (bits == 0) {
      value = left_node(tree, tree[current]).value;
      current = tree[current].subtree_l1[0];
    } else {
      value = right_node(tree, tree[current]).value;
      current = tree[current].subtree_l1[1];
    }
    if (value != NodeValue::NonExisting) {
      // yield data
      return (uint8_t)value;
    }
  }
}

std::vector<int> HuffmanDecoder::HuffmanDecodeDFA(int64_t dc_index, int64_t ac_index, int64_t count_to_read) {
}

