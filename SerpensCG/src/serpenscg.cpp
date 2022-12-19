#include <ap_int.h>
#include <cstdio>
#include <cstring>
#include <cassert>

#include "ap_utils.h"
#include <tapa.h>

#include "serpenscg.h"

constexpr int FIFO_DEPTH = 2;

const int NUM_CH_SPARSE_div_8 = NUM_CH_SPARSE / 8;
const int NUM_CH_SPARSE_mult_4 = NUM_CH_SPARSE * 4;
const int WINDOW_SIZE_div_8 = WINDOW_SIZE / 8;

using double_v4 = tapa::vec_t<double, 4>;

struct MultXVec {
    tapa::vec_t<ap_uint<18>, 4> row;
    double_v4 axv;
};

struct InstRdWr {
    bool rd;
    bool wr;
    //bool require_response;
    int base_addr;
    int len;
};

struct InstVCtrl {
    bool rd;
    bool wr;
    int base_addr;
    int len;
    ap_uint<3> q_rd_idx;
    //ap_uint<3> q_wr_idx;
};

struct InstCmp {
    int M;
    double alpha;
    //ap_uint<3> q_idx;
};

struct ResTerm {
    double res;
    bool term;
};

template <typename T, typename R>
inline void async_read(tapa::async_mmap<T> & A,
                       tapa::ostream<T> & fifo_A,
                       const R i_end_addr,
                       R & i_req,
                       R & i_resp) {
#pragma HLS inline
    if ((i_req < i_end_addr) &
        !A.read_addr.full()) {
        A.read_addr.try_write(i_req);
        ++i_req;
    }
    if (!fifo_A.full() & !A.read_data.empty()) {
        T tmp;
        A.read_data.try_read(tmp);
        fifo_A.try_write(tmp);
        ++i_resp;
    }
}


template <typename T, typename R>
inline void async_write(tapa::async_mmap<T> & Y_out,
                        tapa::istream<T> & fifo_Y,
                        const R num_ite_Y,
                        R & i_req,
                        R & i_resp
                        ) {
#pragma HLS inline
    if ((i_req < num_ite_Y) &
        !fifo_Y.empty() &
        !Y_out.write_addr.full() &
        !Y_out.write_data.full() ) {
        Y_out.write_addr.try_write(i_req);
        T tmpv;
        fifo_Y.try_read(tmpv);
        Y_out.write_data.try_write(tmpv);
        ++i_req;
    }
    uint8_t n_resp;
    if (Y_out.write_resp.try_read(n_resp)) {
        i_resp += R(n_resp) + 1;
    }
}

void rdwr_vec(tapa::async_mmap<double_v8> & vec_p,
              tapa::istream<InstRdWr> & q_inst,
              tapa::istream<double_v8> & q_din,
              tapa::ostream<double_v8> & q_dout,
              tapa::ostream<bool> & q_response
              ) {
    for (;;) {
        auto inst = q_inst.read();
        
        const int rd_end_addr = inst.rd? (inst.base_addr + inst.len) : 0;
        const int wr_end_addr = inst.wr? (inst.base_addr + inst.len) : 0;
        
        const int rd_total = inst.rd? inst.len : 0;
        const int wr_total = inst.wr? inst.len : 0;
        
    rdwr:
        for (int rd_req = inst.base_addr, rd_resp = 0,
             wr_req = inst.base_addr, wr_resp = 0;
             (rd_resp < rd_total) | (wr_resp < wr_total);) {
#pragma HLS loop_tripcount min=1 max=500000
#pragma HLS pipeline II=1
            // rd
            if ((rd_req < rd_end_addr) &
                !vec_p.read_addr.full()) {
                vec_p.read_addr.try_write(rd_req);
                ++rd_req;
            }
            if (!q_dout.full() & !vec_p.read_data.empty()) {
                double_v8 tmp;
                vec_p.read_data.try_read(tmp);
                q_dout.try_write(tmp);
                ++rd_resp;
            }
            
            //wr
            if ((wr_req < wr_end_addr) &
                !q_din.empty() &
                !vec_p.write_addr.full() &
                !vec_p.write_data.full() ) {
                vec_p.write_addr.try_write(wr_req);
                double_v8 tmpv;
                q_din.try_read(tmpv);
                vec_p.write_data.try_write(tmpv);
                ++wr_req;
            }
            uint8_t n_resp;
            if (vec_p.write_resp.try_read(n_resp)) {
                wr_resp += int(n_resp) + 1;
            }
        }
        
        ap_wait();
        
        if (inst.wr){
            q_response.write(true);
        }
    }
}

template <int N>
inline void ctrl_vec(tapa::istream<double_v8> & qm_din,
                     tapa::ostream<double_v8> & qm_dout,
                     tapa::istream<InstVCtrl> & q_inst,
                     tapa::ostream<InstRdWr> & q_mem_inst,
                     tapa::ostreams<double_v8, N> & q_to_pe,
                     tapa::istream<double_v8> & q_updated
                     ) {
#pragma HLS inline
    for(;;) {
#pragma HLS loop_flatten off
        const auto inst_ctrl = q_inst.read();
        
        InstRdWr inst_mem;
        inst_mem.base_addr = inst_ctrl.base_addr;
        inst_mem.len = inst_ctrl.len;
        inst_mem.rd = inst_ctrl.rd;
        inst_mem.wr = inst_ctrl.wr;
        
        const int rd_total = inst_ctrl.rd? inst_ctrl.len : 0;
        const int wr_total = inst_ctrl.wr? inst_ctrl.len : 0;
        
        const int ch_idx = inst_ctrl.q_rd_idx;
        
        q_mem_inst.write(inst_mem);
        ap_wait();
        
    qq:
        for(int i = 0, j = 0; (i < rd_total) | (j < wr_total);) {
#pragma HLS pipeline II=1
            if (!q_updated.empty() & !qm_dout.full()) {
                double_v8 tmp;
                q_updated.try_read(tmp);
                qm_dout.try_write(tmp);
                ++j;
            }
            if (!qm_din.empty() & !q_to_pe[ch_idx].full()) {
                double_v8 tmp;
                qm_din.try_read(tmp);
                q_to_pe[ch_idx].try_write(tmp);
                ++i;
            }
        }
    }
}

inline void ctrl_vec(tapa::istream<double_v8> & qm_din,
                     tapa::ostream<double_v8> & qm_dout,
                     tapa::istream<InstVCtrl> & q_inst,
                     tapa::ostream<InstRdWr> & q_mem_inst,
                     tapa::ostream<double_v8> & q_to_pe,
                     tapa::istream<double_v8> & q_updated
                     ) {
#pragma HLS inline
    for(;;) {
#pragma HLS loop_flatten off
        const auto inst_ctrl = q_inst.read();
        
        InstRdWr inst_mem;
        inst_mem.base_addr = inst_ctrl.base_addr;
        inst_mem.len = inst_ctrl.len;
        inst_mem.rd = inst_ctrl.rd;
        inst_mem.wr = inst_ctrl.wr;
        
        const int rd_total = inst_ctrl.rd? inst_ctrl.len : 0;
        const int wr_total = inst_ctrl.wr? inst_ctrl.len : 0;
        
        q_mem_inst.write(inst_mem);
        ap_wait();
        
    qq:
        for(int i = 0, j = 0; (i < rd_total) | (j < wr_total);) {
#pragma HLS pipeline II=1
            if (!q_updated.empty() & !qm_dout.full()) {
                double_v8 tmp;
                q_updated.try_read(tmp);
                qm_dout.try_write(tmp);
                ++j;
            }
            if (!qm_din.empty() & !q_to_pe.full()) {
                double_v8 tmp;
                qm_din.try_read(tmp);
                q_to_pe.try_write(tmp);
                ++i;
            }
        }
    }
}

void read_edge_list_ptr(tapa::istream<int> & q_inst_in,
                        tapa::ostream<int> & PE_inst,
                        tapa::ostream<int> & rdA_inst,
                        tapa::ostream<int> & aby_inst,
                        tapa::async_mmap<int> & edge_list_ptr
                        ) {
    for (;;) {
#pragma HLS loop_flatten off
        const int num_ite = q_inst_in.read();
        const int M = q_inst_in.read();
        const int A_len = q_inst_in.read();
        
        PE_inst.write(num_ite);
        PE_inst.write(M);
        
        rdA_inst.write(A_len);
        aby_inst.write(M);
        
        const int num_ite_plus1 = num_ite + 1;
        
    rd_ptr:
        for (int i_req = 0, i_resp = 0; i_resp < num_ite_plus1;) {
#pragma HLS loop_tripcount min=1 max=800
#pragma HLS pipeline II=1
            async_read(edge_list_ptr,
                       PE_inst,
                       num_ite_plus1,
                       i_req, i_resp);
        }
    }
}

void read_A(tapa::istream<int> & q_inst_in,
            tapa::ostream<int> & q_inst_out,
            
            tapa::async_mmap<ap_uint<512>> & A,
            tapa::ostream<ap_uint<384>> & fifo_A
            ) {
    for(;;) {
#pragma HLS loop_flatten off
        const int A_len = q_inst_in.read();
        q_inst_out.write(A_len);
        
    rd_A:
        for(int i_req = 0, i_resp = 0; i_resp < A_len;) {
#pragma HLS loop_tripcount min=1 max=10000
#pragma HLS pipeline II=1
            if ((i_req < A_len) &
                !A.read_addr.full()) {
                A.read_addr.try_write(i_req);
                ++i_req;
            }
            if (!fifo_A.full() & !A.read_data.empty()) {
                ap_uint<512> tmp;
                A.read_data.try_read(tmp);
                ap_uint<384> tmp_o;
                for (int i = 0; i < 4; ++i) {
                    tmp_o(i * 96 + 95, i * 96) = tmp(i * 128 + 95, i * 128);
                }
                fifo_A.try_write(tmp_o);
                ++i_resp;
            }
        }
    }
}

void PEG_Xvec(tapa::istream<int> & fifo_inst_in,
              tapa::istream<ap_uint<384>> & fifo_A,
              tapa::istream<double_v8> & fifo_X_in,
              tapa::ostream<int> & fifo_inst_out,
              tapa::ostream<double_v8> & fifo_X_out,
              // to PEG_Yvec
              tapa::ostream<int> & fifo_inst_out_to_Yvec,
              tapa::ostream<MultXVec> & fifo_aXvec
              ) {
    for(;;) {
#pragma HLS loop_flatten off
        const int NUM_ITE = fifo_inst_in.read();
        const int M = fifo_inst_in.read();
        
        fifo_inst_out.write(NUM_ITE);
        fifo_inst_out.write(M);
        
        fifo_inst_out_to_Yvec.write(NUM_ITE);
        fifo_inst_out_to_Yvec.write(M);
        
        double local_X[2][WINDOW_SIZE];
#pragma HLS bind_storage variable=local_X latency=1
#pragma HLS array_partition variable=local_X complete dim=1
#pragma HLS array_partition variable=local_X cyclic factor=X_PARTITION_FACTOR dim=2
        
        auto start_32 = fifo_inst_in.read();
        fifo_inst_out.write(start_32);
        fifo_inst_out_to_Yvec.write(start_32);
        
    main:
        for (int i = 0; i < NUM_ITE; ++i) {
#pragma HLS loop_tripcount min=1 max=49
            
            // fill onchip X
        read_X:
            for (int j = 0; (j < WINDOW_SIZE_div_8) & (j < ((M + 7) >> 3) - i * WINDOW_SIZE_div_8); ) {
#pragma HLS loop_tripcount min=1 max=512
#pragma HLS pipeline II = 1
                if (!fifo_X_in.empty() & !fifo_X_out.full()) {
                    double_v8 x; fifo_X_in.try_read(x);
                    fifo_X_out.try_write(x);
                    for (int kk = 0; kk < 8; ++kk) { //512 / 64 = 8
                        for (int l = 0; l < 2; ++l) {
                            local_X[l][(j << 3) + kk] = x[kk]; // 8 -> << 3
                        }
                    }
                    ++j;
                }
            }
            
            // computation
            const auto end_32 = fifo_inst_in.read();
            fifo_inst_out.write(end_32);
            fifo_inst_out_to_Yvec.write(end_32);
            
        computation:
            for (int j = start_32; j < end_32; ) {
#pragma HLS loop_tripcount min=1 max=200
#pragma HLS pipeline II=1
                if (!fifo_A.empty()) {
                    ap_uint<384> a_pes; fifo_A.try_read(a_pes);
                    MultXVec raxv;
                    
                    for (int p = 0; p < 4; ++p) {
                        ap_uint<96> a = a_pes(95 + p * 96, p * 96);
                        ap_uint<14> a_col = a(95, 82);
                        ap_uint<18> a_row = a(81, 64);
                        ap_uint<64> a_val = a(63,  0);
                        
                        raxv.row[p] = a_row;
                        if (a_row[17] == 0) {
                            double a_val_f64 = tapa::bit_cast<double>(a_val);
                            raxv.axv[p] = a_val_f64 * local_X[p/2][a_col];
                        }
                    }
                    fifo_aXvec.write(raxv);
                    ++j;
                }
            }
            start_32 = end_32;
        }
    }
}

void PEG_Yvec(tapa::istream<int> & fifo_inst_in,
              tapa::istream<MultXVec> & fifo_aXvec,
              tapa::ostream<double> & fifo_Y_out
              ) {
    for(;;) {
#pragma HLS loop_flatten off
        const int NUM_ITE = fifo_inst_in.read();
        const int M = fifo_inst_in.read();
        
        const int num_v_init = (M + NUM_CH_SPARSE_mult_4 - 1) / NUM_CH_SPARSE_mult_4;
        const int num_v_out = (M + NUM_CH_SPARSE - 1) / NUM_CH_SPARSE;
        
        double local_C[4][URAM_DEPTH];
#pragma HLS bind_storage variable=local_C type=RAM_2P impl=URAM latency=1
#pragma HLS array_partition complete variable=local_C dim=1
        
        //init local C
    init_C:
        for (int i = 0; i < num_v_init; ++i) {
#pragma HLS loop_tripcount min=1 max=800
#pragma HLS pipeline II=1
            for (int p = 0; p < 4; ++p) {
                local_C[p][i] = 0.0;
            }
        }
        
        auto start_32 = fifo_inst_in.read();
        
    main:
        for (int i = 0; i < NUM_ITE; ++i) {
#pragma HLS loop_tripcount min=1 max=49
            
            // computation
            const auto end_32 = fifo_inst_in.read();
            
        computation:
            for (int j = start_32; j < end_32; ) {
#pragma HLS loop_tripcount min=1 max=200
#pragma HLS pipeline II=1
#pragma HLS dependence true variable=local_C distance=DEP_DIST_LOAD_STORE
                if (!fifo_aXvec.empty()) {
                    MultXVec raxv; fifo_aXvec.try_read(raxv);
                    for (int p = 0; p < 4; ++p) {
                        auto a_row = raxv.row[p];
                        if (a_row[17] == 0) {
                            local_C[p][a_row] += raxv.axv[p];
                        }
                    }
                    ++j;
                }
            }
            start_32 = end_32;
        }
        
    write_C_outer:
        for (int i = 0, c_idx = 0; i < num_v_out; ++i) {
#pragma HLS loop_tripcount min=1 max=1800
#pragma HLS pipeline II=1
            double out_v = local_C[c_idx][i>>2];
            fifo_Y_out.write(out_v);
            ++c_idx;
            if (c_idx == 4) {c_idx = 0;}
        }
    }
}

void Arbiter_Y(tapa::istream<int> & aby_inst,
               tapa::ostream<int> & aby_inst_out,
               tapa::istreams<double, NUM_CH_SPARSE_div_8> & fifo_in, // 2 = 16 / 8
               tapa::ostream<double> & fifo_out
               ) {
    for(;;) {
#pragma HLS loop_flatten off
        const int M = aby_inst.read();
        aby_inst_out.write(M);
        
        const int num_pe_output = ((M + NUM_CH_SPARSE - 1) / NUM_CH_SPARSE) * NUM_CH_SPARSE_div_8;
        const int num_out = (M + 7) >> 3;
        
    aby:
        for (int i = 0, c_idx = 0; i < num_pe_output;) {
    #pragma HLS loop_tripcount min=1 max=1800
    #pragma HLS pipeline II=1
            if (!fifo_in[c_idx].empty() & !fifo_out.full()) {
                double tmp; fifo_in[c_idx].try_read(tmp);
                if (i < num_out) {
                    fifo_out.try_write(tmp);
                }
                ++i;
                c_idx++;
                if (c_idx == NUM_CH_SPARSE_div_8) {c_idx = 0;}
            }
        }
    }
}

void Merger_Y(tapa::istreams<double, 8> & fifo_in,
              tapa::ostream<double_v8> & fifo_out) {
    for (;;) {
#pragma HLS pipeline II=1
        bool flag_nop = fifo_out.full();
        for (int i = 0; i < 8; ++i) {
            flag_nop |= fifo_in[i].empty();
        }
        if (!flag_nop) {
            double_v8 tmpv;
#pragma HLS aggregate variable=tmpv
            for (int i = 0; i < 8; ++i) {
                double tmp; fifo_in[i].try_read(tmp);
                tmpv[i] = tmp;
            }
            fifo_out.try_write(tmpv);
        }
    }
}

template <typename data_t>
inline void bh(tapa::istream<data_t> & q) {
#pragma HLS inline
    for (;;) {
#pragma HLS pipeline II=1
        data_t tmp; q.try_read(tmp);
    }
}

void black_hole_int(tapa::istream<int> & fifo_in) {
    bh(fifo_in);
}

void black_hole_double_v8(tapa::istream<double_v8> & fifo_in) {
    bh(fifo_in);
}

void black_hole_bool(tapa::istream<bool> & fifo_in) {
    bh(fifo_in);
}

void ctrl_P(tapa::istream<double_v8> & qm_din,
            tapa::ostream<double_v8> & qm_dout,
            tapa::istream<InstVCtrl> & q_inst,
            tapa::ostream<InstRdWr> & q_mem_inst,
            tapa::ostreams<double_v8, 4> & q_to_pe,
            tapa::istream<double_v8> & q_updated
            ) {
    ctrl_vec<4>(qm_din,
                qm_dout,
                q_inst,
                q_mem_inst,
                q_to_pe,
                q_updated);
}

void ctrl_AP(tapa::istream<double_v8> & qm_din,
             tapa::ostream<double_v8> & qm_dout,
             tapa::istream<InstVCtrl> & q_inst,
             tapa::ostream<InstRdWr> & q_mem_inst,
             tapa::ostreams<double_v8, 2> & q_to_pe,
             tapa::istream<double_v8> & q_updated
             ) {
    ctrl_vec<2>(qm_din,
                qm_dout,
                q_inst,
                q_mem_inst,
                q_to_pe,
                q_updated);
}

void ctrl_Z(tapa::istream<double_v8> & qm_din,
            tapa::ostream<double_v8> & qm_dout,
            tapa::istream<InstVCtrl> & q_inst,
            tapa::ostream<InstRdWr> & q_mem_inst,
            tapa::ostreams<double_v8, 2> & q_to_pe,
            tapa::istream<double_v8> & q_updated
            ) {
    ctrl_vec<2>(qm_din,
                qm_dout,
                q_inst,
                q_mem_inst,
                q_to_pe,
                q_updated);
}

void ctrl_X(tapa::istream<double_v8> & qm_din,
            tapa::ostream<double_v8> & qm_dout,
            tapa::istream<InstVCtrl> & q_inst,
            tapa::ostream<InstRdWr> & q_mem_inst,
            tapa::ostream<double_v8> & q_to_pe,
            tapa::istream<double_v8> & q_updated
            ) {
    ctrl_vec(qm_din,
             qm_dout,
             q_inst,
             q_mem_inst,
             q_to_pe,
             q_updated);
}

void ctrl_R(tapa::istream<double_v8> & qm_din,
            tapa::ostream<double_v8> & qm_dout,
            tapa::istream<InstVCtrl> & q_inst,
            tapa::ostream<InstRdWr> & q_mem_inst,
            tapa::ostreams<double_v8, 4> & q_to_pe,
            tapa::istream<double_v8> & q_updated
            ) {
   ctrl_vec<4>(qm_din,
               qm_dout,
               q_inst,
               q_mem_inst,
               q_to_pe,
               q_updated);
}

void read_digA(tapa::istream<int> & q_inst,
               tapa::async_mmap<double_v8> & vec_mem,
               tapa::ostream<double_v8> & q_dout
               ) {
    for(;;) {
#pragma HLS loop_flatten off
        //const int M = q_inst.read();
        const int num_ite = q_inst.read(); // + 7) >> 3;
        
    rd:
        for (int addr_req = 0, i_resp = 0; i_resp < num_ite;) {
#pragma HLS loop_tripcount min=1 max=500000
#pragma HLS pipeline II=1
            async_read(vec_mem,
                       q_dout,
                       num_ite,
                       addr_req, i_resp);
        }
    }
}


/*  computation modules  */

//M2: alpha = rzold / (p' * Ap)
void dot_alpha(tapa::istream<InstCmp> & q_inst,
               tapa::istream<double_v8> & q1,
               tapa::istream<double_v8> & q2,
               tapa::ostream<double> & q_alpha
               ) {
    for(;;) {
#pragma HLS loop_flatten off
        const auto inst = q_inst.read();
        const int M = inst.M;
        const int num_ite = (M + 7) >> 3;
        
        const double rzold = inst.alpha;
        
        double psum[8][DEP_DIST_LOAD_STORE];
#pragma HLS array_partition complete variable=psum dim=1
        
    init:
        for (int i = 0; i < DEP_DIST_LOAD_STORE; ++i) {
#pragma HLS pipeline II=1
            for (int p = 0; p < 8; ++p) {psum[p][i] = 0.0;}
        }
        
    comp1:
        for (int i = 0, idx = 0; i < num_ite; ) {
#pragma HLS pipeline II=1
#pragma HLS dependence true variable=psum distance=DEP_DIST_LOAD_STORE
            if (!q1.empty() & !q2.empty()) {
                double_v8 v1; q1.try_read(v1);
                double_v8 v2; q2.try_read(v2);
                for (int p = 0; p < 8; ++p) {
                    psum[p][idx] += ((i * 8 + p < M)? v1[p] * v2[p] : 0.0);
                }
                ++i;
                ++idx;
                if (idx == DEP_DIST_LOAD_STORE) {idx = 0;}
            }
        }
        
    comp2:
        for (int i = DEP_DIST_LOAD_STORE; i < DEP_DIST_LOAD_STORE * 8; ++i) {
#pragma HLS pipeline II=1
#pragma HLS dependence true variable=psum distance=DEP_DIST_LOAD_STORE
            psum[0][i % DEP_DIST_LOAD_STORE] += psum[i / DEP_DIST_LOAD_STORE][i % DEP_DIST_LOAD_STORE];
        }
        
    comp3:
        for (int i = 1; i < DEP_DIST_LOAD_STORE; ++i) {
#pragma HLS pipeline II=1
            psum[0][0] += psum[0][i];
        }
        
        double pAp = psum[0][0];
        
        double alpha = rzold / pAp;
        
        q_alpha.write(alpha);
    }
}

//M end: res = r' * r
void dot_res(tapa::istream<int> & q_inst,
             tapa::istream<double_v8> & q1,
             tapa::ostream<double> & q2
             ) {
    for(;;) {
#pragma HLS loop_flatten off
        const int M = q_inst.read();
        const int num_ite = (M + 7) >> 3;
        
        double psum[8][DEP_DIST_LOAD_STORE];
#pragma HLS array_partition complete variable=psum dim=1
        
    init:
        for (int i = 0; i < DEP_DIST_LOAD_STORE; ++i) {
#pragma HLS pipeline II=1
            for (int p = 0; p < 8; ++p) {psum[p][i] = 0.0;}
        }
        
    comp1:
        for (int i = 0, idx = 0; i < num_ite; ) {
#pragma HLS pipeline II=1
#pragma HLS dependence true variable=psum distance=DEP_DIST_LOAD_STORE
            if (!q1.empty()) {
                double_v8 v1; q1.try_read(v1);
                for (int p = 0; p < 8; ++p) {
                    psum[p][idx] += ((i * 8 + p < M)? v1[p] * v1[p] : 0.0);
                }
                ++i;
                ++idx;
                if (idx == DEP_DIST_LOAD_STORE) {idx = 0;}
            }
        }
        
    comp2:
        for (int i = DEP_DIST_LOAD_STORE; i < DEP_DIST_LOAD_STORE * 8; ++i) {
#pragma HLS pipeline II=1
#pragma HLS dependence true variable=psum distance=DEP_DIST_LOAD_STORE
            psum[0][i % DEP_DIST_LOAD_STORE] += psum[i / DEP_DIST_LOAD_STORE][i % DEP_DIST_LOAD_STORE];
        }
        
    comp3:
        for (int i = 1; i < DEP_DIST_LOAD_STORE; ++i) {
#pragma HLS pipeline II=1
            psum[0][0] += psum[0][i];
        }
        
        double res = psum[0][0];
        
        q2.write(res);
    }
}

//M6: rznew = r' * z
void dot_rznew(tapa::istream<int> & q_inst,
               tapa::istream<double_v8> & qr,
               tapa::istream<double_v8> & qz,
               tapa::ostream<double> & qrz
               ) {
    for(;;) {
#pragma HLS loop_flatten off
        const int M = q_inst.read();
        const int num_ite = (M + 7) >> 3;
        
        double psum[8][DEP_DIST_LOAD_STORE];
#pragma HLS array_partition complete variable=psum dim=1
        
    init:
        for (int i = 0; i < DEP_DIST_LOAD_STORE; ++i) {
#pragma HLS pipeline II=1
            for (int p = 0; p < 8; ++p) {psum[p][i] = 0.0;}
        }
        
    comp1:
        for (int i = 0, idx = 0; i < num_ite; ) {
#pragma HLS pipeline II=1
#pragma HLS dependence true variable=psum distance=DEP_DIST_LOAD_STORE
            if (!qr.empty() & !qz.empty()) {
                double_v8 v1; qr.try_read(v1);
                double_v8 v2; qz.try_read(v2);
                for (int p = 0; p < 8; ++p) {
                    psum[p][idx] += ((i * 8 + p < M)? v1[p] * v2[p] : 0.0);
                }
                ++i;
                ++idx;
                if (idx == DEP_DIST_LOAD_STORE) {idx = 0;}
            }
        }
        
    comp2:
        for (int i = DEP_DIST_LOAD_STORE; i < DEP_DIST_LOAD_STORE * 8; ++i) {
#pragma HLS pipeline II=1
#pragma HLS dependence true variable=psum distance=DEP_DIST_LOAD_STORE
            psum[0][i % DEP_DIST_LOAD_STORE] += psum[i / DEP_DIST_LOAD_STORE][i % DEP_DIST_LOAD_STORE];
        }
        
    comp3:
        for (int i = 1; i < DEP_DIST_LOAD_STORE; ++i) {
#pragma HLS pipeline II=1
            psum[0][0] += psum[0][i];
        }
        
        double rz = psum[0][0];
        
        qrz.write(rz);
    }
}

//M3: x = x + alpha * p
void updt_x(tapa::istream<InstCmp> & q_inst,
            tapa::istream<double_v8> & qx,
            tapa::istream<double_v8> & qp,
            tapa::ostream<double_v8> & qout
            ) {
    for(;;) {
#pragma HLS loop_flatten off
        const auto inst = q_inst.read();
        const int num_ite = (inst.M + 7) >> 3;

        const double alpha = inst.alpha;

        //qout = x + alpha .* p;
    cc:
        for (int i = 0; i < num_ite;) {
    #pragma HLS pipeline II=1
            if (!qx.empty() & !qp.empty()) {
                double_v8 tmpx; qx.try_read(tmpx);
                double_v8 tmpp; qp.try_read(tmpp);
                qout.write(tmpx + alpha * tmpp);
                ++i;
            }
        }
    }
}

//M7: p = z + (rznew/rzold) * p
void updt_p(tapa::istream<InstCmp> & q_inst,
            
            tapa::istream<double_v8> & qz,
            tapa::istream<double_v8> & qp,
            tapa::ostream<double_v8> & qout
            ) {
    for(;;) {
#pragma HLS loop_flatten off
        const auto inst = q_inst.read();
        const int num_ite = (inst.M + 7) >> 3;

        const double rzndo = inst.alpha;
        
    cc:
        for (int i = 0; i < num_ite;) {
    #pragma HLS pipeline II=1
            if (!qz.empty() & !qp.empty()) {
                double_v8 tmpz; qz.try_read(tmpz);
                double_v8 tmpp; qp.try_read(tmpp);
                qout.write(tmpz + rzndo * tmpp);
                ++i;
            }
        }
    }
}

//M4: r = r - alpha * Ap
void updt_r(tapa::istream<InstCmp> & q_inst,
            tapa::istream<double_v8> & qr,
            tapa::istream<double_v8> & qap,
            tapa::ostream<double_v8> & qout
            ) {
    for(;;) {
#pragma HLS loop_flatten off
        const auto inst = q_inst.read();
        const int num_ite = (inst.M + 7) >> 3;

        const double alpha = inst.alpha;
        
        //qout = x + alpha .* p;
    cc:
        for (int i = 0; i < num_ite;) {
#pragma HLS pipeline II=1
            if (!qr.empty() & !qap.empty()) {
                double_v8 tmpr;   qr.try_read(tmpr);
                double_v8 tmpap; qap.try_read(tmpap);
                qout.write(tmpr - alpha * tmpap);
                ++i;
            }
        }
    }
}

//M5: z = diagA \ r
void left_div(tapa::istream<int> & q_inst,
              tapa::istream<double_v8> & qr,
              tapa::istream<double_v8> & qdiagA,
              tapa::ostream<double_v8> & qz
              ) {
    for(;;) {
#pragma HLS loop_flatten off
        const int M = q_inst.read();
        const int num_ite = (M + 7) >> 3;
    cc:
        for (int i = 0; i < num_ite; ) {
#pragma HLS pipeline II=1
            if (!qr.empty() & !qdiagA.empty()) {
                double_v8 v1; qr.try_read(v1);
                double_v8 v2; qdiagA.try_read(v2);
                double_v8 result;
                for (int p = 0; p < 8; ++p) {
                    result[p] = ((i * 8 + p < M)? v1[p] / v2[p] : 0.0);
                }
                ++i;
                qz.write(result);
            }
        }
    }
}

void wr_r(tapa::async_mmap<double> & vec_r,
          tapa::istream<ResTerm> & q_din
          ) {
    int wr_count = 1;
wr:
    for (int addr_req = 0, i_resp = 0; i_resp < wr_count;) {
#pragma HLS loop_tripcount min=1 max=500000
#pragma HLS pipeline II=1
        if ((addr_req < wr_count) &
            !q_din.empty() &
            !vec_r.write_addr.full() &
            !vec_r.write_data.full() ) {
            vec_r.write_addr.try_write(addr_req);
            ResTerm tmpv;
            q_din.try_read(tmpv);
            vec_r.write_data.try_write(tmpv.res);
            ++addr_req;
            if (tmpv.term) {
                wr_count = addr_req;
            } else {
                ++wr_count;
            }
        }
        uint8_t n_resp;
        if (vec_r.write_resp.try_read(n_resp)) {
            i_resp += int(n_resp) + 1;
        }
    }
}

void global_controller(const int NUM_ITE,
                       const int NUM_A_LEN,
                       const int M,
                       const int rp_time,
                       const double th_termination,
                       
                       tapa::istream<double> & fifo_rr,
                       tapa::istream<double> & fifo_alpha,
                       tapa::istream<double> & fifo_rz,
                       
                       tapa::istream<bool> & resp_P,
                       tapa::istream<bool> & resp_R,
                       tapa::istream<bool> & resp_X,
                       tapa::istream<bool> & resp_Z,
                       tapa::istream<bool> & resp_AP,
                       
                       tapa::ostream<InstVCtrl> & ctrlv_P,
                       tapa::ostream<InstVCtrl> & ctrlv_R,
                       tapa::ostream<InstVCtrl> & ctrlv_X,
                       tapa::ostream<InstVCtrl> & ctrlv_Z,
                       tapa::ostream<InstVCtrl> & ctrlv_AP,
                       tapa::ostream<int> & ctrlv_dA,
                       
                       tapa::ostream<int> & spmv_inst,
                       tapa::ostream<InstCmp> & dotalpha_inst,
                       tapa::ostream<InstCmp> & updtx_inst,
                       tapa::ostream<InstCmp> & updtr_inst,
                       tapa::ostream<int> & leftdiv_inst,
                       tapa::ostream<int> & rz_inst,
                       tapa::ostream<InstCmp> & updtp_inst,
                       tapa::ostream<int> & res_inst,
                       
                       tapa::ostream<ResTerm> & fifo_rr_w
                       ) {
    const int ite_v = (M + 7) >> 3;
    
    double alpha = 1.0;
    double rz_old = 0.0;
    
    for (int rp = -1; rp < rp_time; ++rp) {
#pragma HLS pipeline off
        //M1 - spmv: Ap <- A * p
        ctrlv_P.write({true, false, 0, ite_v, 0});
        ctrlv_AP.write({false, true, 0, ite_v, 0});
        spmv_inst.write(NUM_ITE); spmv_inst.write(M); spmv_inst.write(NUM_A_LEN);
        ap_wait();
        resp_AP.read();
        
        if (rp >= 0) {
            //M2: alpha = rzold / (p' * Ap)
            ctrlv_P.write({true, false, 0, ite_v, 1});
            ctrlv_AP.write({true, false, 0, ite_v, 0});
            dotalpha_inst.write({M, rz_old});
            ap_wait();
            alpha = fifo_alpha.read();
            
            //M3: x = x + alpha * p
            ctrlv_X.write({true, true, 0, ite_v, 0});
            ctrlv_P.write({true, false, 0, ite_v, 2});
            updtx_inst.write({M, alpha});
            ap_wait();
            resp_X.read();
        }
        
        //M4: r = r - alpha * Ap
        ctrlv_R.write({true, true, 0, ite_v, 0});
        ctrlv_AP.write({true, false, 0, ite_v, 1});
        updtr_inst.write({M, alpha});
        ap_wait();
        resp_R.read();
        
        //M end: res = r' * r
        ctrlv_R.write({true, false, 0, ite_v, 3});
        res_inst.write(M);
        ap_wait();
        double res = fifo_rr.read();
        
        bool termination = (res < th_termination) | (rp + 1 == rp_time);
        fifo_rr_w.write({res, termination});
        
        if (termination) break;
        
        //M5: z = diagA \ r
        ctrlv_R.write({true, false, 0, ite_v, 1});
        ctrlv_dA.write(ite_v);
        ctrlv_Z.write({false, true, 0, ite_v, 0});
        leftdiv_inst.write(M);
        ap_wait();
        resp_Z.read();
        
        //M6: rznew = r' * z
        ctrlv_R.write({true, false, 0, ite_v, 2});
        ctrlv_Z.write({true, false, 0, ite_v, 0});
        rz_inst.write(M);
        ap_wait();
        double rz_new = fifo_rz.read();
        
        //M7: p = z + (rznew/rzold) * p
        double rznew_old = (rp < 0)? 0 : (rz_new/rz_old);
        rz_old = rz_new;
        ctrlv_P.write({true, true, 0, ite_v, 3});
        ctrlv_Z.write({true, false, 0, ite_v, 1});
        updtp_inst.write({M, rznew_old});
        ap_wait();
        resp_P.read();
    }
}

void SerpensCG(tapa::mmap<int> edge_list_ptr,
               
               tapa::mmaps<ap_uint<512>, NUM_CH_SPARSE> edge_list_ch,
               
               tapa::mmap<double_v8> vec_x,
               
               tapa::mmap<double_v8> vec_p,
               
               tapa::mmap<double_v8> vec_Ap,
               
               tapa::mmap<double_v8> vec_r,
               
               tapa::mmap<double_v8> vec_z,
               
               tapa::mmap<double_v8> vec_digA,
               
               tapa::mmap<double> vec_res,
               
               const int NUM_ITE,
               const int NUM_A_LEN,
               const int M,
               const int rp_time,
               const double th_termination
               ) {
    //inst fifos
    tapa::stream<int, FIFO_DEPTH> spmv_inst("spmv_inst");
    
    tapa::stream<InstCmp, FIFO_DEPTH> dotalpha_inst("dotalpha_inst");
    
    tapa::stream<InstCmp, FIFO_DEPTH> updtx_inst("updtx_inst");
    
    tapa::stream<InstCmp, FIFO_DEPTH> updtr_inst("updtr_inst");
    
    tapa::stream<int, FIFO_DEPTH> leftdiv_inst("leftdiv_inst");
    
    tapa::stream<int, FIFO_DEPTH> rz_inst("rz_inst");
    
    tapa::stream<InstCmp, FIFO_DEPTH> updtp_inst("updtp_inst");
    
    tapa::stream<int, FIFO_DEPTH> res_inst("res_inst");
    
    // fifos for spmv
    tapa::streams<int, NUM_CH_SPARSE + 1, FIFO_DEPTH> PE_inst("PE_inst");
    
    tapa::streams<int, NUM_CH_SPARSE + 1, FIFO_DEPTH> rdA_inst("rdA_inst");
    
    tapa::streams<int, 8 + 1, FIFO_DEPTH> aby_inst("aby_inst");
    
    tapa::streams<ap_uint<384>, NUM_CH_SPARSE, FIFO_DEPTH> fifo_A("fifo_A");
    
    tapa::streams<int, NUM_CH_SPARSE, FIFO_DEPTH> Yvec_inst("Yvec_inst");
    
    tapa::streams<MultXVec, NUM_CH_SPARSE, FIFO_DEPTH> fifo_aXvec("fifo_aXvec");
    
    tapa::streams<double, NUM_CH_SPARSE, FIFO_DEPTH> fifo_Y_pe("fifo_Y_pe");
    
    tapa::streams<double, 8, FIFO_DEPTH> fifo_Y_pe_abd("fifo_Y_pe_abd");
    
    // fifos for vector modules
    
    //P
    tapa::stream<InstVCtrl, FIFO_DEPTH> fifo_ci_P("fifo_ci_P");
    
    tapa::stream<InstRdWr, FIFO_DEPTH> fifo_mi_P("fifo_mi_P");
    
    tapa::stream<double_v8, FIFO_DEPTH> fifo_din_P("fifo_din_P");
    
    tapa::stream<double_v8, FIFO_DEPTH> fifo_dout_P("fifo_dout_P");
    
    tapa::stream<bool, FIFO_DEPTH> fifo_resp_P("fifo_resp_P");
    
    tapa::streams<double_v8, 4 + NUM_CH_SPARSE, FIFO_DEPTH> fifo_P("fifo_P");
    
    tapa::stream<double_v8, FIFO_DEPTH> fifo_P_updated("fifo_P_updated");
    
    //R
    tapa::stream<InstVCtrl, FIFO_DEPTH> fifo_ci_R("fifo_ci_R");
    
    tapa::stream<InstRdWr, FIFO_DEPTH> fifo_mi_R("fifo_mi_R");
    
    tapa::stream<double_v8, FIFO_DEPTH> fifo_din_R("fifo_din_R");
    
    tapa::stream<double_v8, FIFO_DEPTH> fifo_dout_R("fifo_dout_R");
    
    tapa::stream<bool, FIFO_DEPTH> fifo_resp_R("fifo_resp_R");
    
    tapa::streams<double_v8, 4, FIFO_DEPTH> fifo_R("fifo_R");
    
    tapa::stream<double_v8, FIFO_DEPTH> fifo_R_updated("fifo_R_updated");
    
    //X
    tapa::stream<InstVCtrl, FIFO_DEPTH> fifo_ci_X("fifo_ci_X");
    
    tapa::stream<InstRdWr, FIFO_DEPTH> fifo_mi_X("fifo_mi_X");
    
    tapa::stream<double_v8, FIFO_DEPTH> fifo_din_X("fifo_din_X");
    
    tapa::stream<double_v8, FIFO_DEPTH> fifo_dout_X("fifo_dout_X");
    
    tapa::stream<bool, FIFO_DEPTH> fifo_resp_X("fifo_resp_X");
    
    tapa::stream<double_v8, FIFO_DEPTH> fifo_X("fifo_X");
    
    tapa::stream<double_v8, FIFO_DEPTH> fifo_X_updated("fifo_X_updated");
    
    //AP
    tapa::stream<InstVCtrl, FIFO_DEPTH> fifo_ci_AP("fifo_ci_AP");
    
    tapa::stream<InstRdWr, FIFO_DEPTH> fifo_mi_AP("fifo_mi_AP");
    
    tapa::stream<double_v8, FIFO_DEPTH> fifo_din_AP("fifo_din_AP");
    
    tapa::stream<double_v8, FIFO_DEPTH> fifo_dout_AP("fifo_dout_AP");
    
    tapa::stream<bool, FIFO_DEPTH> fifo_resp_AP("fifo_resp_AP");

    tapa::streams<double_v8, 2, FIFO_DEPTH> fifo_AP("fifo_AP");
    
    tapa::stream<double_v8, FIFO_DEPTH> fifo_AP_updated("fifo_AP_updated");
    
    //Z
    tapa::stream<InstVCtrl, FIFO_DEPTH> fifo_ci_Z("fifo_ci_Z");
    
    tapa::stream<InstRdWr, FIFO_DEPTH> fifo_mi_Z("fifo_mi_Z");
    
    tapa::stream<double_v8, FIFO_DEPTH> fifo_din_Z("fifo_din_Z");
    
    tapa::stream<double_v8, FIFO_DEPTH> fifo_dout_Z("fifo_dout_Z");
    
    tapa::stream<bool, FIFO_DEPTH> fifo_resp_Z("fifo_resp_Z");

    tapa::streams<double_v8, 2, FIFO_DEPTH> fifo_Z("fifo_Z");
    
    tapa::stream<double_v8, FIFO_DEPTH> fifo_Z_updated("fifo_Z_updated");
    
    //diagA
    tapa::stream<int, FIFO_DEPTH> fifo_ci_dA("fifo_ci_dA");
    
    tapa::stream<double_v8, FIFO_DEPTH> fifo_dA("fifo_dA");
    
    //rr
    tapa::stream<double, FIFO_DEPTH> fifo_rr("fifo_rr");
    
    tapa::stream<ResTerm, FIFO_DEPTH> fifo_rr_w("fifo_rr_w");
    
    //alpha
    tapa::stream<double, FIFO_DEPTH> fifo_alpha("fifo_alpha");
    
    //rz
    tapa::stream<double, FIFO_DEPTH> fifo_rz("fifo_rz");

    /* =========deploy modules======= */
    
    tapa::task()
        .invoke<tapa::join>(global_controller,
                            NUM_ITE,
                            NUM_A_LEN,
                            M,
                            rp_time,
                            th_termination,
                            
                            fifo_rr,
                            fifo_alpha,
                            fifo_rz,
                            
                            fifo_resp_P,
                            fifo_resp_R,
                            fifo_resp_X,
                            fifo_resp_Z,
                            fifo_resp_AP,
                            
                            fifo_ci_P,
                            fifo_ci_R,
                            fifo_ci_X,
                            fifo_ci_Z,
                            fifo_ci_AP,
                            fifo_ci_dA,
                            
                            spmv_inst,
                            dotalpha_inst,
                            updtx_inst,
                            updtr_inst,
                            leftdiv_inst,
                            rz_inst,
                            updtp_inst,
                            res_inst,
                            
                            fifo_rr_w
                            )
    
    //P
        .invoke<tapa::detach>(rdwr_vec,
                              vec_p,
                              fifo_mi_P,
                              fifo_dout_P,
                              fifo_din_P,
                              fifo_resp_P
                              )
    
        .invoke<tapa::detach>(ctrl_P,
                              fifo_din_P,
                              fifo_dout_P,
                              fifo_ci_P,
                              fifo_mi_P,
                              fifo_P,
                              fifo_P_updated
                              )
    
    //R
        .invoke<tapa::detach>(rdwr_vec,
                              vec_r,
                              fifo_mi_R,
                              fifo_dout_R,
                              fifo_din_R,
                              fifo_resp_R
                              )
    
        .invoke<tapa::detach>(ctrl_R,
                              fifo_din_R,
                              fifo_dout_R,
                              fifo_ci_R,
                              fifo_mi_R,
                              fifo_R,
                              fifo_R_updated
                              )
    
    //X
        .invoke<tapa::detach>(rdwr_vec,
                              vec_x,
                              fifo_mi_X,
                              fifo_dout_X,
                              fifo_din_X,
                              fifo_resp_X
                              )
    
        .invoke<tapa::detach>(ctrl_X,
                              fifo_din_X,
                              fifo_dout_X,
                              fifo_ci_X,
                              fifo_mi_X,
                              fifo_X,
                              fifo_X_updated
                              )
    
    //AP
        .invoke<tapa::detach>(rdwr_vec,
                              vec_Ap,
                              fifo_mi_AP,
                              fifo_dout_AP,
                              fifo_din_AP,
                              fifo_resp_AP
                              )
    
        .invoke<tapa::detach>(ctrl_AP,
                              fifo_din_AP,
                              fifo_dout_AP,
                              fifo_ci_AP,
                              fifo_mi_AP,
                              fifo_AP,
                              fifo_AP_updated
                              )
    
    //Z
        .invoke<tapa::detach>(rdwr_vec,
                              vec_z,
                              fifo_mi_Z,
                              fifo_dout_Z,
                              fifo_din_Z,
                              fifo_resp_Z
                              )
    
        .invoke<tapa::detach>(ctrl_Z,
                              fifo_din_Z,
                              fifo_dout_Z,
                              fifo_ci_Z,
                              fifo_mi_Z,
                              fifo_Z,
                              fifo_Z_updated
                              )
    
    //digA
        .invoke<tapa::detach>(read_digA,
                              fifo_ci_dA,
                              vec_digA,
                              fifo_dA
                              )
    
    //SpMV - 1
        .invoke<tapa::detach>(read_edge_list_ptr,
                              spmv_inst,
                              PE_inst,
                              rdA_inst,
                              aby_inst,
                              edge_list_ptr
                              )
    
        .invoke<tapa::detach, NUM_CH_SPARSE>(read_A,
                                             rdA_inst,
                                             rdA_inst,
                                             edge_list_ch,
                                             fifo_A
                                             )
        .invoke<tapa::detach>(black_hole_int,
                              rdA_inst)
    
        .invoke<tapa::detach>(PEG_Xvec,
                              PE_inst,
                              fifo_A,
                              fifo_P,
                              PE_inst,
                              fifo_P,
                              Yvec_inst,
                              fifo_aXvec
                              )
                        
        .invoke<tapa::detach, NUM_CH_SPARSE>(PEG_Yvec,
                                             Yvec_inst,
                                             fifo_aXvec,
                                             fifo_Y_pe
                                             )
    
        .invoke<tapa::detach, 8>(Arbiter_Y,
                                 aby_inst,
                                 aby_inst,
                                 fifo_Y_pe,
                                 fifo_Y_pe_abd
                                 )
        .invoke<tapa::detach>(black_hole_int,
                              aby_inst)
    
        .invoke<tapa::detach>(Merger_Y,
                              fifo_Y_pe_abd,
                              fifo_AP_updated
                              )
    
    //M2: alpha = rzold / (p' * Ap)
        .invoke<tapa::detach>(dot_alpha,
                              dotalpha_inst,
                              fifo_P,
                              fifo_AP,
                              fifo_alpha
                              )
    
    //M3: x = x + alpha * p
        .invoke<tapa::detach>(updt_x,
                              updtx_inst,
                              fifo_X,
                              fifo_P,
                              fifo_X_updated
                              )
    
    //M4: r = r - alpha * Ap
        .invoke<tapa::detach>(updt_r,
                              updtr_inst,
                              fifo_R,
                              fifo_AP,
                              fifo_R_updated
                              )
    
    //M5: z = diagA \ r
        .invoke<tapa::detach>(left_div,
                              leftdiv_inst,
                              fifo_R,
                              fifo_dA,
                              fifo_Z_updated
                              )
    
    //M6: rznew = r' * z
        .invoke<tapa::detach>(dot_rznew,
                              rz_inst,
                              fifo_R,
                              fifo_Z,
                              fifo_rz
                              )
    
    //M7: p = z + (rznew/rzold) * p
        .invoke<tapa::detach>(updt_p,
                              updtp_inst,
                              fifo_Z,
                              fifo_P,
                              fifo_P_updated
                              )
    
    //M residual
        .invoke<tapa::detach>(dot_res,
                              res_inst,
                              fifo_R,
                              fifo_rr
                              )
    
        .invoke<tapa::join>(wr_r,
                            vec_res,
                            fifo_rr_w
                            )
    
    //SpMV - 2
        .invoke<tapa::detach, NUM_CH_SPARSE - 1>(PEG_Xvec,
                                                 PE_inst,
                                                 fifo_A,
                                                 fifo_P,
                                                 PE_inst,
                                                 fifo_P,
                                                 Yvec_inst,
                                                 fifo_aXvec
                                                 )
        .invoke<tapa::detach>(black_hole_int,
                              PE_inst)
        .invoke<tapa::detach>(black_hole_double_v8,
                              fifo_P)
    ;
}
