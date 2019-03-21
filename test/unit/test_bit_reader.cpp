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
  BitReader bs{"\xd7"};
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
  BitReader bs{str};
  ASSERT_EQ(bs.read_bit(), 1);
  ASSERT_EQ(bs.read_bit(), 1);
  ASSERT_EQ(bs.read_bit(), 1);
  ASSERT_EQ(bs.read_bit(), 1);
}

TEST_F(TestBitReader, TestReadBitsFF) {
  string str{"\xff\00", 2};
  BitReader bs{str};
  int bits;
  bs.read_nbits(7, bits);
  ASSERT_EQ(bits, 0b01111111);
}

TEST_F(TestBitReader, TestReadBits) {
  string str{"\xcf\00", 2};
  BitReader bs{str};
  int bits;
  bs.read_nbits(7, bits);
  ASSERT_EQ(bits, 0b1100111);
}

TEST_F(TestBitReader, TestReturnBit) {
  string str{"\xcf\01", 2};
  BitReader bs{str};
  int bits;
  bs.read_nbits(7, bits);
  ASSERT_EQ(bits, 0b1100111);
  bs.read_nbits(2, bits);
  ASSERT_EQ(bits, 0b10);

  bs.return_bit();
  bs.return_bit();
  bs.read_nbits(2, bits);
  ASSERT_EQ(bits, 0b10);
}

