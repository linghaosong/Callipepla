#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <cassert>

#include <ap_int.h>
#include <tapa.h>

#include "../include/mmio.h"
#include "../include/sparse_helper.h"
#include "../SerpensCG/src/serpenscg.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::vector;
using std::min;
using std::max;

template <typename T>
using aligned_vector = std::vector<T, tapa::aligned_allocator<T>>;

int main(int argc, char **argv) {
    gflags::ParseCommandLineFlags(&argc, &argv, /*remove_flags=*/true);
    
    printf("start host JPCG\n");
    int num_ites = 1;

    if (argc == 3) {
        num_ites = atoi(argv[2]);
    }
    else if (argc != 2) {
        cout << "Usage: " << argv[0] << " [matrix A file] [num_ites]" << std::endl;
        return EXIT_FAILURE;
    }

    char * filename_A = argv[1];

    cout << "num_ites = " << num_ites <<  "\n";
    
    FILE * fin = fopen(filename_A, "rb");
    
    int M;
    fread(&M, sizeof(int), 1, fin);
    
    int sparse_A_fpga_vec_size;
    fread(&sparse_A_fpga_vec_size, sizeof(int), 1, fin);
    
    int MAX_SIZE_edge_LIST_PTR;
    fread(&MAX_SIZE_edge_LIST_PTR, sizeof(int), 1, fin);
    
    int MAX_LEN_edge_PTR;
    fread(&MAX_LEN_edge_PTR, sizeof(int), 1, fin);
    
    int edge_list_ptr_fpga_chunk_size;
    fread(&edge_list_ptr_fpga_chunk_size, sizeof(int), 1, fin);
    
    
    
    int vec_fpga_column_size = ((M + 8 - 1) / 8) * 8;
    int vec_fpga_chunk_size = ((vec_fpga_column_size + 511)/512) * 512;
    
    
    aligned_vector<double> diag_A_fpga(vec_fpga_chunk_size, 0.0);
    for (int j = 0; j < M; ++j) {
        double tmpx;
        fread(&tmpx, sizeof(double), 1, fin);
        diag_A_fpga[j] = tmpx;
    }
    
    aligned_vector<int> edge_list_ptr_fpga(edge_list_ptr_fpga_chunk_size, 0);
    for (int j = 0; j < edge_list_ptr_fpga_chunk_size; ++j) {
        int tmpx;
        fread(&tmpx, sizeof(int), 1, fin);
        edge_list_ptr_fpga[j] = tmpx;
    }
    
    vector<aligned_vector<ap_uint<128>> > sparse_A_fpga_vec(NUM_CH_SPARSE);
    for (int i = 0; i < NUM_CH_SPARSE; ++i) {
        sparse_A_fpga_vec[i].resize(sparse_A_fpga_vec_size, 0);
        for (int j = 0; j < sparse_A_fpga_vec_size; ++j) {
            ap_uint<128> tmp128;
            fread(&tmp128, sizeof(ap_uint<128>), 1, fin);
            sparse_A_fpga_vec[i][j] = tmp128;
        }
    }
    
    fclose(fin);
    
    /* CG */
    vector<double> x(M, 0.0);
    vector<double> b(M, 1.0);

    vector<double> r(M, 0.0);
    vector<double> p(M, 0.0);
    vector<double> z(M, 0.0);
    
    //copy x to p to perform the preprocessing step of CG solver
    for (int i = 0; i < M; ++i) {
        p[i] = x[i];
    }
    
    //copy b to r to perform the preprocessing step of CG solver
    for (int i = 0; i < M; ++i) {
        r[i] = b[i];
    }
    
    aligned_vector<double> vec_X_fpga(vec_fpga_chunk_size, 0.0);
    aligned_vector<double> vec_P_fpga(vec_fpga_chunk_size, 0.0);
    aligned_vector<double> vec_AP_fpga(vec_fpga_chunk_size, 0.0);
    aligned_vector<double> vec_R_fpga(vec_fpga_chunk_size, 0.0);
    aligned_vector<double> vec_Z_fpga(vec_fpga_chunk_size, 0.0);
    
    for (int i = 0; i < M; ++i) {
        vec_X_fpga[i] = x[i];
        vec_P_fpga[i] = p[i];
        vec_R_fpga[i] = r[i];
    }
    
    aligned_vector<double> vec_RES_fpga(((num_ites + 1 + 511)/512) * 512, 0.0);
    
    std::string bitstream;
    if (const auto bitstream_ptr = getenv("TAPAB")) {
        bitstream = bitstream_ptr;
    }
    
    const double th_termination = 1e-12;
    
    cout << "launch kernel\n";
    double time_taken
    = tapa::invoke(SerpensCG, bitstream,
                   tapa::read_only_mmap<int>(edge_list_ptr_fpga),
                   tapa::read_only_mmaps<ap_uint<128>, NUM_CH_SPARSE>(sparse_A_fpga_vec).reinterpret<ap_uint<512> >(),
                   
                   tapa::read_write_mmap<double>(vec_X_fpga).reinterpret<double_v8>(),
                   tapa::read_write_mmap<double>(vec_P_fpga).reinterpret<double_v8>(),
                   tapa::read_write_mmap<double>(vec_AP_fpga).reinterpret<double_v8>(),
                   tapa::read_write_mmap<double>(vec_R_fpga).reinterpret<double_v8>(),
                   tapa::read_write_mmap<double>(vec_Z_fpga).reinterpret<double_v8>(),
                   
                   tapa::read_only_mmap<double>(diag_A_fpga).reinterpret<double_v8>(),
                   
                   tapa::write_only_mmap<double>(vec_RES_fpga),
                   
                   MAX_SIZE_edge_LIST_PTR,
                   MAX_LEN_edge_PTR,
                   M,
                   num_ites,
                   th_termination
                   );
    int ite_kernel;
    for (ite_kernel = num_ites;
         (ite_kernel > 0) && (vec_RES_fpga[ite_kernel] < 1e-305);
         --ite_kernel);
    
    cout << "res form device ..." << endl;
    for (int i = 0; i < ite_kernel + 1; ++i) {
        cout << "Ite = " << std::setw(5) << i << ", rr = " << vec_RES_fpga[i] << endl;
    }
    
    time_taken *= 1e-9; // total time in second
    printf("Kernel time is %f ms\n", time_taken*1000);
    
    cout << "# Iterations = " << ite_kernel << endl;
    
    time_taken /= ite_kernel;
    printf("Per iteration time is %f ms\n", time_taken*1000);
    
    return 0;
}
