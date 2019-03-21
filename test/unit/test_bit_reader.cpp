#include <gtest/gtest.h>
#include <bit_reader.hpp>
#include <string>
#include <sstream>
#include <iostream>

using namespace std;
class TestBitReader :public ::testing::Test {
public:



};

TEST_F(TestBitReader, TestReadOneBit) {
  // 0xff will be escaped by a following 0x00
  string str{"\xd7"};
  std::stringstream ss{str};
  BitReader bs{ss};
  ASSERT_EQ(bs.read_bit(), 1);
  ASSERT_EQ(bs.read_bit(), 1);
  ASSERT_EQ(bs.read_bit(), 0);
  ASSERT_EQ(bs.read_bit(), 1);
  ASSERT_EQ(bs.read_bit(), 0);
  ASSERT_EQ(bs.read_bit(), 1);
  ASSERT_EQ(bs.read_bit(), 1);
  ASSERT_EQ(bs.read_bit(), 1);
}

TEST_F(TestBitReader, TestReadOneBitFF) {
  // 0xff will be escaped by a following 0x00
  string str{"\xff\00", 2};
  std::stringstream ss{str};
  BitReader bs{ss};
  ASSERT_EQ(bs.read_bit(), 1);
  ASSERT_EQ(bs.read_bit(), 1);
  ASSERT_EQ(bs.read_bit(), 1);
  ASSERT_EQ(bs.read_bit(), 1);
}

TEST_F(TestBitReader, TestReadBitsFF) {
  string str{"\xff\00", 2};
  std::stringstream ss{str};
  BitReader bs{ss};
  ASSERT_EQ(bs.read_nbits(7), 0b01111111);
}

TEST_F(TestBitReader, TestReadBits) {
  string str{"\xcf\00", 2};
  std::stringstream ss{str};
  BitReader bs{ss};
  ASSERT_EQ(bs.read_nbits(7), 0b1100111);
}

