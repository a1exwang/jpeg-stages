#include <vector>
#include "huffman.hpp"
#include <string>
#include <iostream>
#include "decode.hpp"
#include <glog/logging.h>

using namespace std;


constexpr int64_t MaximumHuffmanDepth = 16;

int extend(int value, int nBits) {
  if (value >> (nBits - 1) == 0) {
    value -= (1 << nBits) - 1;
  }
  return value;
}

std::vector<int>
HuffmanDecoder::HuffmanDecode(int64_t dc_index, int64_t ac_index, int64_t count_to_read) {
  auto read_nbits = [this](int n) {
    int ret = 0;
    for (int i = 0; i < n; i++) {
      ret |= (read_bit() << (n-i-1));
    }
    return ret;
  };

  auto read_tree = [this](std::vector<int64_t> &tree) -> uint8_t {
    int current = 0;
    while (true) {
      int bit = read_bit();
      if (bit < 0) {
        cerr << "Invalid entropy stream, unexpected EOF" << endl;
        abort();
      }
      // 0 for left, 1 for right
      if (bit == 0) {
        // left child = 2*k + 1
        current = 2 * current+1;
      } else {
        // right child = 2*k + 2
        current = 2 * current+2;
      }
      if (current >= tree.size()) {
        cout << "Invalid entropy stream, code too long" << endl;
        abort();
      }
      if (tree[current] != -1) {
        // yield data
        return (uint8_t)tree[current];
      }
    }
  };

  vector<int> ret(count_to_read, 0);
  int64_t count_read = 0;
  // Finish building huffman tree
  auto t = read_tree(dc_trees[dc_index]);
  int diff = 0;
  if (t == 0) {
    diff = 0;
  } else {
    diff = read_nbits(t);
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
    count_read += r;

    // decode_zz k
    int zzk = read_nbits(ssss);
    ret[count_read++] = extend(zzk, ssss);
  }
  return ret;
}

int HuffmanDecoder::read_bit() {
  if (remaining_bits == 0) {
    octet = is.get();
    if (octet == 0xFF) {
      auto b = is.get();
      CHECK(b == 0);
    }
    if (!is) {
      return -1;
    }
    remaining_bits = 8;
  }
  // Take the MSB(Most significant bit)
  auto ret = (octet & (1 << (remaining_bits-1))) ? 1 : 0;
  remaining_bits--;
  return ret;
}

HuffmanDecoder::HuffmanDecoder(
    std::istream &is,
    std::map<int64_t, HuffmanTable> dc_dhts,
    std::map<int64_t, HuffmanTable> ac_dhts
) :is(is), dc_trees(dc_dhts.size(), std::vector<int64_t>(1<<17, -1)), ac_trees(ac_dhts.size(), std::vector<int64_t>(1<<17, -1)) {


  // Build Huffman trees
  for (int i = 0; i < dc_dhts.size(); i++) {
    auto &dc_values = dc_dhts[i].codes;
    auto &dc_counts = dc_dhts[i].counts;
    auto &dc_tree = dc_trees[i];
    auto it = dc_values.begin();
    int64_t skip = 0;
    for (int depth = 0; depth < dc_counts.size(); depth++) {
      skip = (skip + (depth == 0 ? 0 : dc_counts[depth-1]))* 2;
      for (int j = 0; j < /*width*/ dc_counts[depth]; j++) {
        dc_tree[(1 << (depth+1)) - 1 + j + skip] = *it++;
      }
    }
  }

  for (int i = 0; i < ac_dhts.size(); i++) {
    auto &ac_values = ac_dhts[i].codes;
    auto &ac_counts = ac_dhts[i].counts;
    auto &ac_tree = ac_trees[i];
    auto it = ac_values.begin();
    int skip = 0;
    for (int depth = 0; depth < ac_counts.size(); depth++) {
      skip = (skip + (depth == 0 ? 0 : ac_counts[depth-1]))* 2;
      for (int j = 0; j < /*width*/ ac_counts[depth]; j++) {
        ac_tree[(1 << (depth+1)) - 1 + j + skip] = *it++;
      }
    }
  }
}

