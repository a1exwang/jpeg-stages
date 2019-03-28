#include "decode.hpp"
#include "huffman.hpp"
#include <iostream>
#include <fstream>
#include <atomic>
#include <thread>
#include <numeric>
#include <streambuf>
#include <string>
#include <iterator>
#include <map>
#include <algorithm>
#include <fftw3.h>
#include <cmath>
#include <glog/logging.h>

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>

#include <opencv2/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <immintrin.h>
#include <cstdlib>

using namespace std;

constexpr int kBlockSide = 8;
constexpr int kBlockSize = kBlockSide * kBlockSide;
using IntTable = std::array<int, kBlockSize>;
using FloatTable = std::array<float, kBlockSize>;
constexpr uint8_t JpegSectionStartByte = 0xff;
typedef enum {kDc = 0, kAc = 1, kNumCoefficientKinds = 2} CoefficientKind;

template <typename T>
T clip(T x, T lb, T ub) {
  return std::min(std::max(x, lb), ub);
}
constexpr float fac = 0.5f / kBlockSide;
constexpr float fac0 = 1.4142135623730951f * fac;
constexpr float fac00 = 1.0f / kBlockSide;


constexpr uint8_t JpegSOF0Mark = 0xc0;
constexpr uint8_t JpegDHTMark = 0xc4;
constexpr uint8_t JpegSOIMark = 0xd8;
constexpr uint8_t JpegEOIMark = 0xd9;
constexpr uint8_t JpegSOSMark = 0xda;
constexpr uint8_t JpegDQTMark = 0xdb;
constexpr uint8_t JpegAPP0Mark = 0xe0;
constexpr uint8_t JpegCommentMark = 0xfe;

int zigzag_to_normal_table[] = {
    0,  1,  8, 16,  9,  2,  3, 10,
    17, 24, 32, 25, 18, 11,  4,  5,
    12, 19, 26, 33, 40, 48, 41, 34,
    27, 20, 13,  6,  7, 14, 21, 28,
    35, 42, 49, 56, 57, 50, 43, 36,
    29, 22, 15, 23, 30, 37, 44, 51,
    58, 59, 52, 45, 38, 31, 39, 46,
    53, 60, 61, 54, 47, 55, 62, 63
};

template <class InputIterator, class OutputIterator>
static void fill_matrix_in_zigzag(InputIterator input, OutputIterator output, int nRows, int nCols) {
  int i = 0;
  int j = 0;
  int di = -1;
  int dj = 1;
  for (int diag = 0; diag < nRows + nCols - 1; ++diag) {
    while (i >= 0 && j >= 0 && i < nRows && j < nCols) {
      *(output + i * nCols + j) = *input;
      ++input;
      i += di;
      j += dj;
    }
    i -= di;
    j -= dj;
    di *= -1;
    dj *= -1;
    if ((i == 0 || i == nRows - 1) && j < nCols - 1) {
      ++j;
    } else {
      ++i;
    }
  }
}

template <class InputIterator, class OutputIterator>
static void fill_matrix_in_zigzag_fast64(InputIterator input, OutputIterator output) {
  for (int i = 0; i < 64; i++) {
    output[zigzag_to_normal_table[i]] = input[i];
  }
}

map<uint8_t, string> MarkerNames = {
    {JpegSOF0Mark, "SOF Start of Frame 0 Baseline DCT"},
    {JpegDHTMark, "DHT Define Huffman Table"},
    {JpegSOIMark, "SOI Start of Image"},
    {JpegEOIMark, "EOI End of Image"},
    {JpegSOSMark, "SOS Start Of Scan"},
    {JpegDQTMark, "DQT Define Quantitization Table"},
    {JpegAPP0Mark, "APP0"},
    {JpegCommentMark, "Comment"},
};

struct JFIFApp0Section {
  uint8_t length_high;
  uint8_t length_low;
  char jfif_str[5];
  uint8_t version[2];
  uint8_t density_unit;
  uint16_t x_density;
  uint16_t y_density;
  uint8_t x_thumbnail;
  uint8_t y_thumbnail;
  char thumbnail_data[0]; // 3 * x_thumbnail * y_thumbnail
  bool is_valid() const {
    return jfif_str == string("JFIF");
  }
} __attribute__((packed));

struct StartOfFrame0 {
  uint16_t sample_precision;
  uint16_t lines, cols;
  uint16_t component_count;
  struct ComponentParameter {
    uint16_t component_id;
    uint16_t horizontal_sampling_factor;
    uint16_t vertical_sampling_factor;
    uint16_t quantization_table_id;
  };
  uint16_t hmax, vmax;
  vector<ComponentParameter> component_parameters;
};

class QuantizationTable {
public:
  int8_t precision; // always 0 for Baseline DCT
  int8_t id;
  array<uint16_t, 64> q;
};


class JpegEOF : public std::exception {
public:
  const char *what() const noexcept { return what_.c_str(); }

private:
  string what_ = "EOF";
};

class JpegFileData {
public:
  JpegFileData(string data, int64_t decode_batch_size) :jpeg_data(move(data)), decode_batch_size(decode_batch_size) { }

  void ParseSections() {
    while (true) {
      // start of section
      auto b = read_byte();
      CHECK(b == JpegSectionStartByte);
      uint8_t marker = read_byte();
      switch (marker) {
        case JpegSOIMark: {
          break;
        }
        case JpegAPP0Mark: {
          app0 = *(JFIFApp0Section*)current_data();
          CHECK(app0.is_valid());
          auto thumbnail_data_size = 3 * app0.x_thumbnail * app0.y_thumbnail;
          skip(sizeof(JFIFApp0Section) + thumbnail_data_size);
          break;
        }
        case JpegEOIMark: {
          if (offset != jpeg_data.size()) {
            LOG(WARNING) << "Trailing data after EOI size=" << jpeg_data.size() - offset;
          }
//          LOG(INFO) << "File reading done";
          return;
          break;
        }
        case JpegDQTMark: {
//          LOG(WARNING) << "DQT";
          auto l = read_be<uint16_t>();
          l -= 2;
          while (l > 0) {
            auto q = parse_q();
            quantization_tables.resize(max((size_t)q.id + 1, quantization_tables.size()));
            quantization_tables[q.id] = q;
            // zigzag currently
            l -= 65;
          }
          break;
        }
        case JpegSOSMark: {
          // length
          read_be<uint16_t>();
          auto n_components = read_byte();

          for (int i = 0; i < n_components; i++) {
            auto i_component = read_byte();
            CHECK(i_component == i + 1);
            auto b = read_byte();
            auto dc_id = b >> 4;
            auto ac_id = b & 0xfu;
            ac_ids[i] = ac_id;
            dc_ids[i] = dc_id;
//            LOG(WARNING) << "dc, ac: " << dc_id << ", " << ac_id;
          }
          auto ss = read_byte();
          auto se = read_byte();
          auto ahai = read_byte();
          CHECK(ss == 0);
          CHECK(se == 63);
          CHECK(ahai == 0);

          decoded_data = read_entropy();
          for (int c = 0; c < sof0.component_count; c++) {
            dense_decoded_data.emplace_back(mcu_col_count*mcu_row_count*kBlockSize, 0);
            for (int i = 0; i < mcu_row_count; i++) {
              for (int j = 0; j < mcu_col_count; j++) {
                for (int ii = 0; ii < kBlockSide; ii++) {
                  for (int jj = 0; jj < kBlockSide; jj++) {
                    auto left_index = i*mcu_col_count*kBlockSide*kBlockSide +
                                      j*kBlockSide*kBlockSide +
                                      ii * kBlockSide +
                                      jj;
                    dense_decoded_data[c][left_index] =
                        decoded_data[c][i*mcu_col_count+j][ii * kBlockSide + jj];
                  }
                }
              }
            }
          }
          break;
        }
        case JpegDHTMark: {
          auto dht = parse_dht();
          if (dht.is_ac) {
            ac_dhts[dht.id] = dht;
          } else {
            dc_dhts[dht.id] = dht;
          }
          break;
        }
        case JpegSOF0Mark: {
          sof0 = parse_sof0();
          stringstream ss;
          ss << "SOF0: " << sof0.cols << "x" << sof0.lines << ", P: " << sof0.sample_precision << "bits, ";
          ss << sof0.component_count << " components, ";
          for (const auto &p : sof0.component_parameters) {
            ss << "<"
               << p.horizontal_sampling_factor << ":"
               << p.vertical_sampling_factor << ","
               << p.quantization_table_id << ">";
          }
//          LOG(WARNING) << ss.str();
          break;
        }
        default: {
          auto length = skip_section();
//          LOG(WARNING) << "Skipping section with marker 0x" << hex << (uint64_t)marker
//              << (MarkerNames.find(marker) != MarkerNames.end() ? "("+MarkerNames[marker]+")" : "")
//              << dec << " size " << length;

        }
      }
    }
  }

  // total size is width*height*batch_size
  void apply_idct_batch_inplace(float *table, float *output, int64_t width, int64_t height, int64_t batch_size) {
    CHECK(width == 8 && height == 8);

    // prepare
    int64_t block_size = width * height;

    __m256 first_row_fac = _mm256_set_ps(fac0, fac0, fac0, fac0, fac0, fac0, fac0, fac00);
    __m256 row_fac = _mm256_set_ps(fac, fac, fac, fac, fac, fac, fac, fac0);
    for (int64_t b = 0; b < batch_size; b++) {
      int k = 0;
      int64_t batch_offset = b*block_size;
      // per block
      for (int i = 0; i < width; ++i) {
        __m256 row = _mm256_load_ps(&table[batch_offset + k]);
        if (i == 0) {
          row = _mm256_mul_ps(row, first_row_fac);
        } else {
          row = _mm256_mul_ps(row, row_fac);
        }
        _mm256_store_ps(&table[batch_offset + k], row);
        k += 8;
      }
    }

    int n[] = {kBlockSide, kBlockSide};
    fftwf_r2r_kind kinds[] = {FFTW_REDFT01, FFTW_REDFT01};
    fftwf_plan plan = fftwf_plan_many_r2r(2, n, batch_size, table, n, 1, block_size, output, n, 1, block_size, kinds, 0);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
  }

  void DecodeDataBatch() {
    for (int c = 0; c < sof0.component_count; c++) {

      auto my_decoded_data = std::unique_ptr<float[], decltype(&free)>((float*)aligned_alloc(256, kBlockSize * mcu_count*sizeof(float)), &free);
      const auto &current_q_table =
          quantization_tables[sof0.component_parameters[c].quantization_table_id];
      for (int x = 0; x < mcu_count; x++) {
        for (int i = 0; i < kBlockSize; i++) {
          my_decoded_data[x*kBlockSize + i] = decoded_data[c][x][i] * current_q_table.q[i];
        }
      }
      done_batched.push_back(
          std::unique_ptr<float, decltype(&free)>((float*)aligned_alloc(256, kBlockSize * mcu_count*sizeof(float)), &free)
      );
      apply_idct_batch_inplace(my_decoded_data.get(), done_batched[c].get(), kBlockSide, kBlockSide, mcu_count);
    }
  }

  cv::Mat ConvertColorSpaceAVX();
public:
  QuantizationTable parse_q() {
    uint8_t qh = read_byte();
    CHECK((qh >> 4) == 0);
    auto id = qh & 0xf;
    QuantizationTable ret;
    ret.id = id;
    ret.precision = 0;

    array<uint8_t, ret.q.size()> tmp;
    for (int i = 0; i < ret.q.size(); i++) {
      auto b = read_byte();
      tmp[i] = b;
    }
//    fill_matrix_in_zigzag(tmp.begin(), ret.q.begin(), kBlockSide, kBlockSide);
    fill_matrix_in_zigzag_fast64(tmp.begin(), ret.q.begin());
    return ret;
  }

  StartOfFrame0 parse_sof0() {
    auto length = read_be<uint16_t>();
    StartOfFrame0 ret;
    ret.sample_precision = read_byte();
    ret.lines = read_be<uint16_t>();
    ret.cols = read_be<uint16_t>();
    ret.component_count = read_byte();

    uint16_t hmax = 0, vmax = 0;
    for (int i = 0; i < ret.component_count; i++) {
      StartOfFrame0::ComponentParameter p;
      p.component_id = read_byte();
      uint8_t a = read_byte();
      p.vertical_sampling_factor = a & (uint8_t)0xFU;
      p.horizontal_sampling_factor = a >> (uint8_t)4U;
      p.quantization_table_id = read_byte();
      ret.component_parameters.push_back(p);
      if (p.vertical_sampling_factor > vmax) {
        vmax = p.vertical_sampling_factor;
      }
      if (p.horizontal_sampling_factor > hmax) {
        hmax = p.horizontal_sampling_factor;
      }
    }
    ret.hmax = hmax;
    ret.vmax = vmax;
    return ret;
  };

  HuffmanTable parse_dht() {
    // length
    read_be<uint16_t>();
    HuffmanTable table;
    auto a = read_byte();
    table.is_ac = (a >> 4);
    table.id = a & 0xf;
    stringstream ss;
    int x = 0;
    table.counts.resize(16);
    for (auto &i : table.counts) {
      i = read_byte();
      ss << (size_t)i << " ";
      x += i;
    }
//    LOG(WARNING) << "DHT: " << table.id << ", " << ss.str() << " total: " << x;
    for (int i = 0; i < table.counts.size(); i++) {
      table.table.emplace_back();
      for (int j = 0;j < table.counts[i]; j++) {
        auto c = read_byte();
        table.table[i].push_back(c);
        table.codes.push_back(c);
      }
    }
    return table;
  }
  int64_t skip_section() {
    auto length = read_be <uint16_t >();

    skip(length - 2);
    return length;
  }
  int64_t decode_batch_size;
  std::vector<std::vector<IntTable>> read_entropy() {
    string str1 = jpeg_data.substr(offset);

    std::vector<std::vector<IntTable>> tables(sof0.component_count);

    auto r = kBlockSide * sof0.vmax;
    auto c = kBlockSide * sof0.hmax;
    mcu_row_count = ((sof0.lines + r - 1) / r);
    mcu_col_count = ((sof0.cols + c - 1) / c);
    mcu_count = mcu_row_count * mcu_col_count;

    HuffmanDecoder decoder(move(str1), dc_dhts, ac_dhts, decode_batch_size);

    int64_t total_saved = 0;
    for (int c = 0; c < sof0.component_count; ++c) {
      std::vector<IntTable> &component_tables = tables[c];
      int n_tables =
          sof0.component_parameters[c].horizontal_sampling_factor *
          sof0.component_parameters[c].vertical_sampling_factor;
      component_tables.reserve(n_tables * mcu_count);
    }
    for (int mcu = 0; mcu < mcu_count; ++mcu) {
      for (int c = 0; c < sof0.component_count; ++c) {
        std::vector<IntTable>& component_tables = tables[c];

        int n_tables =
            sof0.component_parameters[c].horizontal_sampling_factor *
            sof0.component_parameters[c].vertical_sampling_factor;
        for (int t = 0; t < n_tables; ++t) {
          IntTable table;

          auto data = decoder.HuffmanDecode(dc_ids[c], ac_ids[c], kBlockSize);
//          fill_matrix_in_zigzag(data.begin(), table.begin(), kBlockSide, kBlockSide);
          fill_matrix_in_zigzag_fast64(data.data(), table.data());

          if (!component_tables.empty()) {
            // differential encoded for DC values
            table[0] += component_tables.back()[0];
          }
          component_tables.push_back(table);
        }
      }
    }
    original_size = jpeg_data.size();
    after_huffman = (3*sof0.cols*sof0.lines - total_saved);
    final_size = (3*sof0.cols*sof0.lines);
//    skip(ss1.tellg());
    skip(decoder.get_offset());
    return tables;
  }
  uint8_t read_byte() {
    if (offset >= jpeg_data.size()) {
      throw JpegEOF();
    }
    return (uint8_t)jpeg_data[offset++];
  }
  template <typename T>
  T read_be() {
    uint8_t tmp[sizeof(T)];
    for (int i = 0; i < sizeof(T); i++) {
      tmp[sizeof(T) - i - 1] = current_data()[i];
    }
    skip(sizeof(T));

    return *reinterpret_cast<const T*>(tmp);
  }
  const uint8_t *current_data() const {
    return (const uint8_t*)&jpeg_data[offset];
  }
  void skip(int64_t n) {
    offset += n;
  }

  string jpeg_data;

  int64_t offset = 0;

  JFIFApp0Section app0;
  StartOfFrame0 sof0;
  map<int64_t, HuffmanTable> ac_dhts;
  map<int64_t, HuffmanTable> dc_dhts;
  map<uint64_t, uint64_t> ac_ids, dc_ids;

  vector<QuantizationTable> quantization_tables;

  vector<vector<float>> dense_decoded_data;
  vector<vector<IntTable>> decoded_data;
  vector<std::unique_ptr<float, decltype(&free)>> done_batched;

  int mcu_count, mcu_row_count, mcu_col_count;
  int64_t original_size;
  int64_t after_huffman;
  int64_t final_size;
};


cv::Mat JpegFileData::ConvertColorSpaceAVX() {
  cv::Mat mat((int)ceil(double(sof0.lines)/kBlockSide)*kBlockSide, (int)ceil(double(sof0.cols)/kBlockSide)*kBlockSide, CV_8UC3);
  for (int i = 0; i < mcu_row_count; i++) {
    for (int j = 0; j < mcu_col_count; j++) {
      for (int x = 0; x < kBlockSide; x++) {
        float *yyptr = &done_batched.at(0).get()[(i * mcu_col_count + j) * kBlockSize + x * kBlockSide];
        float *cbptr = &done_batched.at(1).get()[(i * mcu_col_count + j) * kBlockSize + x * kBlockSide];
        float *crptr = &done_batched.at(2).get()[(i * mcu_col_count + j) * kBlockSize + x * kBlockSide];

        __m256 yyvec = _mm256_load_ps(yyptr);
        __m256 rvec = _mm256_set_ps(128.0f, 128.0f, 128.0f, 128.0f, 128.0f, 128.0f, 128.0f, 128.0f);

        {
          __m256 crcbvec = _mm256_load_ps(crptr);
          __m256 factorvec = _mm256_set_ps(1.402f,1.402f,1.402f,1.402f,1.402f,1.402f,1.402f,1.402f);
          crcbvec = _mm256_mul_ps(crcbvec, factorvec);
          rvec = _mm256_add_ps(rvec, yyvec);
          rvec = _mm256_add_ps(rvec, crcbvec);
        }

        __m256 bvec = _mm256_set_ps(128.0f, 128.0f, 128.0f, 128.0f, 128.0f, 128.0f, 128.0f, 128.0f);
        {
          __m256 cb1vec = _mm256_load_ps(cbptr);
          __m256 cb1factorvec = _mm256_set_ps(1.772f,1.772f,1.772f,1.772f,1.772f,1.772f,1.772f,1.772f);
          cb1vec = _mm256_mul_ps(cb1vec, cb1factorvec);
          bvec = _mm256_add_ps(bvec, yyvec);
          bvec = _mm256_add_ps(bvec, cb1vec);
        }

        __m256 gvec = _mm256_set_ps(128.0f, 128.0f, 128.0f, 128.0f, 128.0f, 128.0f, 128.0f, 128.0f);
        {
          __m256 crcbvec = _mm256_load_ps(crptr);
          __m256 crcbfactorvec = _mm256_set_ps(-0.71414f,-0.71414f,-0.71414f,-0.71414f,-0.71414f,-0.71414f,-0.71414f,-0.71414f);
          crcbvec = _mm256_mul_ps(crcbvec, crcbfactorvec);
          gvec = _mm256_add_ps(gvec, crcbvec);

          crcbvec = _mm256_load_ps(cbptr);
          crcbfactorvec = _mm256_set_ps(-0.34414f,-0.34414f,-0.34414f,-0.34414f,-0.34414f,-0.34414f,-0.34414f,-0.34414f);
          crcbvec = _mm256_mul_ps(crcbvec, crcbfactorvec);
          gvec = _mm256_add_ps(gvec, crcbvec);

          gvec = _mm256_add_ps(gvec, yyvec);
        }

        __m256 zerovec = _mm256_set_ps(0,0,0,0,0,0,0,0);
        __m256 maxvec = _mm256_set_ps(255,255,255,255,255,255,255,255);
        rvec = _mm256_max_ps(rvec, zerovec);
        gvec = _mm256_max_ps(gvec, zerovec);
        bvec = _mm256_max_ps(bvec, zerovec);

        rvec = _mm256_min_ps(rvec, maxvec);
        gvec = _mm256_min_ps(gvec, maxvec);
        bvec = _mm256_min_ps(bvec, maxvec);
        __m256i rveci = _mm256_cvtps_epi32(rvec);
        __m256i gveci = _mm256_cvtps_epi32(gvec);
        __m256i bveci = _mm256_cvtps_epi32(bvec);

        auto row = i * kBlockSide + x;
        auto col = j * kBlockSide;
        for (int y = 0; y < kBlockSide; y++, col++) {
          auto offset = row * mat.cols*3 + col*3;
          ((uint8_t*)mat.data)[offset] = (uint8_t)((int32_t*)&bveci)[y];
          ((uint8_t*)mat.data)[offset+1] = (uint8_t)((int32_t*)&gveci)[y];
          ((uint8_t*)mat.data)[offset+2] = (uint8_t)((int32_t*)&rveci)[y];
        }
      }
    }
  }
  cv::Rect myROI(0, 0, sof0.cols, sof0.lines);
  return mat(myROI);

}

cv::Mat jst_decode(const std::string &jpeg_data, int64_t decode_batch_size, bool verbose) {
  JpegFileData jpeg(jpeg_data, decode_batch_size);

  jpeg.ParseSections();
  jpeg.DecodeDataBatch();

  if (verbose) {
    LOG(WARNING)
        << "Original file size: " << jpeg.original_size/1e3 << "KB, "
        << "After Huffman decoding: " << jpeg.after_huffman/1e3 << "KB,"
        << " (" << 100.0*jpeg.after_huffman/jpeg.original_size << "% of original, " << 100.0*jpeg.after_huffman/jpeg.final_size <<  "% of final), "
        << "Matrix size: " << jpeg.final_size/1e3 << "KB, (" << 100.0*jpeg.final_size/jpeg.original_size << "%)";
  }
  return jpeg.ConvertColorSpaceAVX();
}

