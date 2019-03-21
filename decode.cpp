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

using namespace std;

const int kBlockSide = 8;
const int kBlockSize = kBlockSide * kBlockSide;
using IntTable = std::array<int, kBlockSize>;
using DoubleTable = std::array<double, kBlockSize>;
constexpr uint8_t JpegSectionStartByte = 0xff;
typedef enum {kDc = 0, kAc = 1, kNumCoefficientKinds = 2} CoefficientKind;

template <typename T>
T clip(T x, T lb, T ub) {
  return std::min(std::max(x, lb), ub);
}


constexpr uint8_t JpegSOF0Mark = 0xc0;
constexpr uint8_t JpegDHTMark = 0xc4;
constexpr uint8_t JpegSOIMark = 0xd8;
constexpr uint8_t JpegEOIMark = 0xd9;
constexpr uint8_t JpegSOSMark = 0xda;
constexpr uint8_t JpegDQTMark = 0xdb;
constexpr uint8_t JpegAPP0Mark = 0xe0;
constexpr uint8_t JpegCommentMark = 0xfe;

template <class InputIterator, class OutputIterator>
void fill_matrix_in_zigzag(InputIterator input, OutputIterator output, int nRows, int nCols) {
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
  JpegFileData(string data) :jpeg_data(move(data)) { }

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

  void applyIDCT(const IntTable& table, DoubleTable &output) {
    DoubleTable input;

    int k = 0;
    constexpr double fac = 0.5 / kBlockSide;
    double fac0 = std::sqrt(2);
    for (int i = 0; i < kBlockSide; ++i) {
      for (int j = 0; j < kBlockSide; ++j) {
        input[k] = table[k] * fac;
        if (i == 0) {
          input[k] *= fac0;
        }
        if (j == 0) {
          input[k] *= fac0;
        }
        ++k;
      }
    }

    fftw_plan plan = fftw_plan_r2r_2d(kBlockSide, kBlockSide, input.data(), output.data(),
                                      FFTW_REDFT01, FFTW_REDFT01, 0);

    fftw_execute(plan);
    fftw_destroy_plan(plan);
  }
  // total size is width*height*batch_size
  void apply_idct_batch(const double *table, double *output, int64_t width, int64_t height, int64_t batch_size) {
    auto input = make_unique<double[]>(static_cast<size_t>(width * height * batch_size));

    // prepare
    int64_t block_size = width * height;
    for (int64_t b = 0; b < batch_size; b++) {
      int k = 0;
      int64_t batch_offset = b*block_size;
      constexpr double fac = 0.5 / kBlockSide;
      double fac0 = std::sqrt(2);
      for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
          input[batch_offset + k] = table[batch_offset + k] * fac;
          if (i == 0) {
            input[batch_offset + k] *= fac0;
          }
          if (j == 0) {
            input[batch_offset + k] *= fac0;
          }
          ++k;
        }
      }
    }
//    for (int64_t b = 0; b < batch_size; b++) {
//      int64_t batch_offset = b*block_size;
//      fftw_plan plan = fftw_plan_r2r_2d(kBlockSide, kBlockSide, input.get() + batch_offset, output + batch_offset,
//                                        FFTW_REDFT01, FFTW_REDFT01, 0);
//
//      fftw_execute(plan);
//      fftw_destroy_plan(plan);
//    }
    int n[] = {kBlockSide, kBlockSide};
    fftw_r2r_kind kinds[] = {FFTW_REDFT01, FFTW_REDFT01};
    fftw_plan plan = fftw_plan_many_r2r(2, n, batch_size, input.get(), n, 1, block_size, output, n, 1, block_size, kinds, 0);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
  }

  void DecodeDataBatch() {
    for (int c = 0; c < sof0.component_count; c++) {
      done.emplace_back(mcu_count);

      auto my_decoded_data = make_unique<double[]>(static_cast<size_t>(kBlockSide * kBlockSide * mcu_count));
      const auto &current_q_table =
          quantization_tables[sof0.component_parameters[c].quantization_table_id];
      for (int x = 0; x < mcu_count; x++) {
        for (int i = 0; i < kBlockSize; i++) {
          my_decoded_data[x*kBlockSize + i] = decoded_data[c][x][i] * current_q_table.q[i];
        }
      }
      done_batched.emplace_back(kBlockSize*mcu_count, 0);
      apply_idct_batch(my_decoded_data.get(), done_batched[c].data(), kBlockSide, kBlockSide, mcu_count);
    }
  }

  cv::Mat ConvertColorSpaceBatched() {
    cv::Mat mat(sof0.lines, sof0.cols, CV_8UC3);
    for (int i = 0; i < mcu_row_count; i++) {
      for (int j = 0; j < mcu_col_count; j++) {
        for (int x = 0; x < kBlockSide; x++) {
          for (int y = 0; y < kBlockSide; y++) {
            auto row = i * kBlockSide + x;
            auto col = j * kBlockSide + y;

            double yy = done_batched[0][(i*mcu_col_count + j)*kBlockSize + x * kBlockSide + y];
            double cb = done_batched[1][(i*mcu_col_count + j)*kBlockSize + x * kBlockSide + y];
            double cr = done_batched[2][(i*mcu_col_count + j)*kBlockSize + x * kBlockSide + y];


            double r = yy + 1.402f * cr + 128.0f;
            double g = yy - 0.34414f * cb - 0.71414f * cr + 128.0f;
            double b = yy + 1.772f * cb + 128.0f;

            b = std::floor(clip(b, 0.0, 255.0));
            g = std::floor(clip(g, 0.0, 255.0));
            r = std::floor(clip(r, 0.0, 255.0));
            if (row < sof0.lines && col < sof0.cols) {
              mat.at<cv::Vec3b>(row, col) = cv::Vec3b{(uint8_t)b, (uint8_t)g, (uint8_t)r};
            }
          }
        }
      }
    }
    return mat;
  }
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
    fill_matrix_in_zigzag(tmp.begin(), ret.q.begin(), kBlockSide, kBlockSide);
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
  std::vector<std::vector<IntTable>> read_entropy() {
    string str1 = jpeg_data.substr(offset);

    std::vector<std::vector<IntTable>> tables(sof0.component_count);

    auto r = kBlockSide * sof0.vmax;
    auto c = kBlockSide * sof0.hmax;
    mcu_row_count = ((sof0.lines + r - 1) / r);
    mcu_col_count = ((sof0.cols + c - 1) / c);
    mcu_count = mcu_row_count * mcu_col_count;

    HuffmanDecoder decoder(move(str1), dc_dhts, ac_dhts);

    int64_t total_saved = 0;
    for (int mcu = 0; mcu < mcu_count; ++mcu) {
      for (int c = 0; c < sof0.component_count; ++c) {
        std::vector<IntTable>& component_tables = tables[c];

        int n_tables =
            sof0.component_parameters[c].horizontal_sampling_factor *
            sof0.component_parameters[c].vertical_sampling_factor;
        for (int t = 0; t < n_tables; ++t) {
          IntTable table;

          auto data = decoder.HuffmanDecode(dc_ids[c], ac_ids[c], kBlockSize);
          fill_matrix_in_zigzag(data.begin(), table.begin(), kBlockSide, kBlockSide);

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

  vector<vector<double>> dense_decoded_data;
  vector<vector<IntTable>> decoded_data;
  vector<vector<DoubleTable>> done;
  // c * (mcu_count*kBlockSize)
  vector<vector<double>> done_batched;

  int mcu_count, mcu_row_count, mcu_col_count;
  int64_t original_size;
  int64_t after_huffman;
  int64_t final_size;
};

cv::Mat jst_decode(const std::string &jpeg_data, bool verbose) {
  JpegFileData jpeg(jpeg_data);

  jpeg.ParseSections();
  jpeg.DecodeDataBatch();

  if (verbose) {
    LOG(WARNING)
        << "Original file size: " << jpeg.original_size/1e3 << "KB, "
        << "After Huffman decoding: " << jpeg.after_huffman/1e3 << "KB,"
        << " (" << 100.0*jpeg.after_huffman/jpeg.original_size << "% of original, " << 100.0*jpeg.after_huffman/jpeg.final_size <<  "% of final), "
        << "Matrix size: " << jpeg.final_size/1e3 << "KB, (" << 100.0*jpeg.final_size/jpeg.original_size << "%)";
  }
  return jpeg.ConvertColorSpaceBatched();
}

