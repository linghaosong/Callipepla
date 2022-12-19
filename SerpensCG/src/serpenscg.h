#ifndef SEEPENSCG_H
#define SEEPENSCG_H

#include <ap_int.h>
#include <tapa.h>

constexpr int NUM_CH_SPARSE = 16;

constexpr int X_PARTITION_FACTOR = 4; // BRAMs = 512 / 16 / 2 = 16 -> factor = 16 / (64 / 16)
constexpr int WINDOW_SIZE = X_PARTITION_FACTOR * 1024;
constexpr int DEP_DIST_LOAD_STORE = 7;
constexpr int URAM_DEPTH = 3 * 4096 * 2;
// ch16: 3 * 4096 * 2 * 16 * 4 = 1,572,864

using double_v8 = tapa::vec_t<double, 8>;

void SerpensCG(tapa::mmap<int> edge_list_ptr,
               tapa::mmaps<ap_uint<512>, NUM_CH_SPARSE> edge_list_ch,
               tapa::mmap<double_v8> vec_x,
               tapa::mmap<double_v8> vec_p,
               tapa::mmap<double_v8> vec_Ap,
               tapa::mmap<double_v8> vec_r,
               tapa::mmap<double_v8> vec_z,
               tapa::mmap<double_v8> vec_digA,
               tapa::mmap<double> vec_res,
               const int NUM_ITE, const int NUM_A_LEN, const int M,
               const int rp_time,
               const double th_termination
               );

#endif  // SEEPENSCG_H
