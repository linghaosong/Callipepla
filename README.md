[![DOI](https://zenodo.org/badge/547376034.svg)](https://zenodo.org/badge/latestdoi/547376034)

# Callipepla & SerpensCG

Callipepla & SerpensCG are Conjugate Gradient (CG) solvers on High Bandwdith Memory (HBM) FPGA (Xilinx Alveo U280). 

Requirements: 
+ TAPA + Autobridge
+ Following [Install TAPA](https://tapa.readthedocs.io/en/release/installation.html) to install TAPA(Autobridge) and Gurobi.
+ Vitis 2021.2
+ Xilinx xilinx_u280_xdma_201920_3 shell and a Xilinx U280 FPGA card.

### Input matrix format & sample input
The host code takes martrix market format(https://math.nist.gov/MatrixMarket/formats.html). We test on sparse matrices from SuiteSparse(https://sparse.tamu.edu) collection. We have a few local examples in the [matrices](https://github.com/linghaosong/Callipepla/tree/main/matrices) folder.

## To build host:
For SerpensCG, 

    cd SerpensCG; mkdir build; cd build;
    g++ -o serpenscg -Wno-write-strings -Wunused-result -O2 ../src/serpenscg.cpp ../src/serpenscg-host.cpp -ltapa -lfrt -lglog -lgflags -lOpenCL 

For Callipepla, 

    cd Callipepla; mkdir build; cd build;
    g++ -o callipepla -Wno-write-strings -Wunused-result -O2 ../src/callipepla.cpp ../src/callipepla-host.cpp -ltapa -lfrt -lglog -lgflags -lOpenCL 

## To run on board:
For SerpensCG, 

    cd SerpensCG
    TAPAB=./bitstream/SerpensCG_xilinx_u280_xdma_201920_3.xclbin ./serpenscg ../matrices/mhd3200b/mhd3200b.mtx 100

For Callipepla, 

    cd Callipepla
    TAPAB=./bitstream/Callipepla_xilinx_u280_xdma_201920_3.xclbin ./callipepla ../matrices/mhd3200b/mhd3200b.mtx 100



To learn more about the techinqual details, please see [Callipepla: Stream Centric Instruction Set and Mixed Precision for Accelerating Conjugate Gradient Solver](https://arxiv.org/abs/2209.14350).


If you find this code useful, please cite:

    @inproceedings{song2023callipepla,
        title={Callipepla: Stream Centric Instruction Set and Mixed Precision for Accelerating Conjugate Gradient Solver},
        author={Song, Linghao and Guo, Licheng and Basalama, Suhail and Chi, Yuze and Lucas, Robert F and Cong, Jason},
        booktitle={The 2023 ACM/SIGDA International Symposium on Field-Programmable Gate Arrays},
        year = {2023}
    }
