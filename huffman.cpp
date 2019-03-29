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

std::vector<int> HuffmanDecoder::HuffmanDecode(int64_t dc_index, int64_t ac_index, size_t count_to_read) {
  vector<int> ret(count_to_read, 0);
  size_t count_read = 0;
  // Finish building huffman tree
  bool done;
  int jump;
  int zzk = read_tree(dc_trees[dc_index], true, done, jump);
  ret.at(count_read++) = zzk;

  while (count_read < count_to_read) {
    zzk = read_tree(ac_trees.at(ac_index), false, done, jump);
    if (done) {
      break;
    }
    count_read += jump;
    ret.at(count_read++) = zzk;
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
  Node *current = nullptr;
  while (true) {
    int bit = bit_reader.read_bit();
    if (bit < 0) {
      cerr << "Invalid entropy stream, unexpected EOF" << endl;
      abort();
    }
    // 0 for left, 1 for right

    NodeValue value = current->children[bit]->value;
    current = current->children[bit];
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


    // root node
    int n = 1;
    for (int depth = 0; depth < ac_counts.size(); depth++) {
      skip = (skip + (depth == 0 ? 0 : ac_counts[depth - 1])) * 2;
      // leaf nodes
      n += ac_counts[depth];
      // inner nodes
      n += (1UL << (depth + 1)) - (skip + ac_counts[depth]);
    }
    tree.resize(n);

    // root node
    size_t index = 0;

    tree[index].value = NodeValue::NonExisting;
    tree[index].index = index;
    index++;

    size_t last_mid_node_start = 0;
    skip = 0;
    for (int depth = 0; depth < ac_counts.size(); depth++) {
      skip = (skip + (depth == 0 ? 0 : ac_counts[depth-1]))* 2;
      // leaf nodes
      for (int j = 0; j < /*width*/ ac_counts[depth]; j++) {
        tree.at(index).value = NodeValue(*it++);
        tree.at(index).index = index;
        index++;
        // <depth+1, skip+j>

        // if last_mid_node left is empty
        //   put in left
        // else
        //   put in right
        //   last_mid_node_start++;
        if (tree.at(last_mid_node_start).children[0] == nullptr) {
          tree.at(last_mid_node_start).children[0] = &tree.at(index-1);
        } else {
          tree.at(last_mid_node_start).children[1] = &tree.at(index-1);
          last_mid_node_start++;
        }
      }
      auto mid_node_start = index;
      // inner nodes
      for (size_t j = skip + ac_counts[depth]; j < (1UL << (depth+1)); j++) {
        tree[index].value = NodeValue::NonExisting;
        tree[index].index = index;
        index++;
        if (tree.at(last_mid_node_start).children[0] == nullptr) {
          tree.at(last_mid_node_start).children[0] = &tree.at(index-1);
        } else {
          tree.at(last_mid_node_start).children[1] = &tree.at(index-1);
          last_mid_node_start++;
        }
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

int HuffmanDecoder::read_tree_batched(vector<HuffmanDecoder::Node> &tree, bool is_dc, bool &done, int &jump) {
  int zzk = 0;
  Node *current = &tree[0];
  while (true) {
    int bits;
    int bit_count = bit_reader.read_nbits(bit_batch_size, bits);
    int original_bit_count = bit_count;
//    if (bits < 0) {
//      cerr << "Invalid entropy stream, unexpected EOF" << endl;
//      abort();
//    }
//    CHECK(bit_count > 0);
    uint8_t value;
    bool is_leaf = read_tree_nbits(&current, bit_count, bits, value);
    if (bit_count > 0 && !is_leaf) {
      bit_reader.return_nbits(bit_count);
    }
    if (is_leaf) {
      int ssss = is_dc ? value : (value & 0xF);
      zzk = 0;

      static int bit_count_mask[] = {
          0b0, 0b1, 0b11, 0b111,
          0b1111, 0b11111, 0b111111, 0b1111111,
          0b11111111, 0b111111111, 0b1111111111, 0b11111111111,
          0b111111111111, 0b1111111111111, 0b11111111111111, 0b111111111111111,
      };
      // TODO(aocheng): Determine which implementation is faster
      // Current benchmarks show the table implemention is slightly faster(0.97% vs 0.70% of total time)
      // But that might depends on other facters, so I leave the other implementation here for further tests.
//      bits &= (1<<bit_count) - 1;
      bits &= bit_count_mask[bit_count];
      if (bit_count > ssss) {
        bit_reader.return_nbits(bit_count - ssss);
        zzk = bits >> (bit_count - ssss);
      } else if (bit_count < ssss) {
        bit_reader.read_nbits(ssss - bit_count, zzk);
        zzk |= bits << (ssss - bit_count);
      } else {
        zzk = bits;
      }
      zzk = extend(zzk, ssss);
      jump = (value >> 4) & 0xF;

      done = ssss == 0 && jump != 15;
      return zzk;
    }
  }
}

bool HuffmanDecoder::read_tree_nbits(HuffmanDecoder::Node **current, int &batch_size, int bits, uint8_t &value) {
  auto next_node = (*current)->subtrees[batch_size-1][bits];
  /**
   * Make one big step in the Huffman tree.
   * There are 3 circumstances:
   * 1. The step is too big, let's go half of the step size;
   * 2. The step is too small to reach a leaf node, return current node and bit count;
   * 3. The step is just enough to reach a leaf node, return the leaf value.
   */
  if (next_node == nullptr) {
//    CHECK(batch_size > 1);
    auto shift_size = (batch_size-batch_size/2);
    batch_size /= 2;
    bool ok = read_tree_nbits(current, batch_size, bits >> shift_size, value);
    batch_size += shift_size;
    return ok;
  } else {
    auto val = next_node->value;
    if (int64_t(val) >= 0) {
      batch_size = 0;
      value = int8_t(val);
      return true;
    } else {
      batch_size = 0;
      *current = next_node;
      return false;
    }
  }
}


HuffmanDecoder::Node *HuffmanDecoder::take_precedence(
    vector<HuffmanDecoder::Node> &tree,
    HuffmanDecoder::Node &node,
    int64_t child_id,
    int64_t nbits) {

  CHECK(nbits <= 63 && nbits > 0);
  Node *current_node = &node;
  for (int64_t i = nbits-1; i >= 0; i--) {
    int64_t current_bit = (child_id >> i) & 0x1;
    current_node = current_node->children[current_bit];
//    CHECK(next_index >= 0);
    if (current_node == nullptr) {
      return nullptr;
    }
  }
  return current_node;

}

