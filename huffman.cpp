#include <vector>
#include "huffman.hpp"
#include <string>
#include <iostream>
#include "decode.hpp"
#include <glog/logging.h>

using namespace std;


constexpr int64_t MaximumHuffmanDepth = 16;

static int extend(int value, int nBits) {
  // value < 2^nBits
  if (value >> (nBits - 1) == 0) {
    value -= (1 << nBits) - 1;
  }
  return value;
}

std::vector<int>
HuffmanDecoder::HuffmanDecode(int64_t dc_index, int64_t ac_index, int64_t count_to_read) {
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

  ret[count_read++] = diff;

  while (count_read < count_to_read) {
    auto ac_value = read_tree(ac_trees[ac_index]);
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
    ret[count_read++] = extend(zzk, ssss);
  }
  return ret;
}

HuffmanDecoder::HuffmanDecoder(
    std::istream &is,
    std::map<int64_t, HuffmanTable> dc_dhts,
    std::map<int64_t, HuffmanTable> ac_dhts
) :bit_reader(is), dc_trees(dc_dhts.size(), std::vector<Node>(1<<17, {NodeValue::NonExisting})), ac_trees(ac_dhts.size(), std::vector<Node>(1<<17, {NodeValue::NonExisting})) {

  convert_table_to_tree(dc_dhts, dc_trees);
  convert_table_to_tree(ac_dhts, ac_trees);
//  build_subtree(dc_trees);
//  build_subtree(ac_trees);

}

uint8_t HuffmanDecoder::read_tree(std::vector<Node> &tree) {
  size_t current = 0;
  while (true) {
    int bit = bit_reader.read_bit();
    if (bit < 0) {
      cerr << "Invalid entropy stream, unexpected EOF" << endl;
      abort();
    }
    // 0 for left, 1 for right
    if (bit == 0) {
      // left child = 2*k + 1
      current = (current<<1)+1;
    } else {
      // right child = 2*k + 2
      current = (current<<1)+2;
    }
    if (current >= tree.size()) {
      cout << "Invalid entropy stream, code too long" << endl;
      abort();
    }
    if (tree.at(current).value != NodeValue::NonExisting) {
      // yield data
      return (uint8_t)tree.at(current).value;
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
    int skip = 0;
    for (int depth = 0; depth < ac_counts.size(); depth++) {
      skip = (skip + (depth == 0 ? 0 : ac_counts[depth-1]))* 2;
      for (int j = 0; j < /*width*/ ac_counts[depth]; j++) {
        Node node{NodeValue(*it++)};
        tree[(1 << (depth+1)) - 1 + j + skip] = node;
      }
    }
  }

}

void HuffmanDecoder::build_subtree(vector<vector<HuffmanDecoder::Node>> &trees) {
  for (auto &tree : trees) {
    for (int i = 0; i < tree.size(); i++) {
      if (tree[i].value != NodeValue::NonExisting) {
        if (2*i+2 < tree.size()) {
          tree[i].subtree_l1[0] = tree[2*i+1].value;
          tree[i].subtree_l1[1] = tree[2*i+2].value;
        }

        if (2*(2*i+2)+2 < tree.size()) {
          tree[i].subtree_l2[0b00] = tree[2*(2*i+1)+1].value;
          tree[i].subtree_l2[0b01] = tree[2*(2*i+1)+2].value;
          tree[i].subtree_l2[0b10] = tree[2*(2*i+2)+1].value;
          tree[i].subtree_l2[0b11] = tree[2*(2*i+2)+2].value;
        }
      }
    }
  }

}

uint8_t HuffmanDecoder::read_tree_batched(vector<HuffmanDecoder::Node> &tree) {
  auto step = [&tree](int current, int bit) {
    // 0 for left, 1 for right
    if (bit == 0) {
      // left child = 2*k + 1
      current = (current<<1)+1;
    } else {
      // right child = 2*k + 2
      current = (current<<1)+2;
    }
    if (current >= tree.size()) {
      cout << "Invalid entropy stream, code too long" << endl;
      abort();
    }
    return current;
  };
  size_t current = 0;
  while (true) {
    int bits;
    int nread = bit_reader.read_nbits(2, bits);
    if (nread <= 0) {
      cerr << "Invalid entropy stream, unexpected EOF" << endl;
      abort();
    }
    if (nread == 2) {
      auto subtree_value = tree.at(current).subtree_l1[bits];
      if (int64_t(subtree_value) >= 0) {
        return uint8_t(int64_t(subtree_value));
      } else {

      }
    } else {
      // step nread
      CHECK(nread == 2);
    }

    if (int64_t(tree.at(current).value) >= 0) {
      // yield data
      return (uint8_t)tree.at(current).value;
    }
  }
}

std::vector<int> HuffmanDecoder::HuffmanDecodeDFA(int64_t dc_index, int64_t ac_index, int64_t count_to_read) {
}

