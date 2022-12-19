# Callipepla & SerpensCG

Callipepla & SerpensCG are Conjugate Gradient (CG) solvers on High Bandwdith Memory (HBM) FPGA (Xilinx Alveo U280). 

Requirements: 
+ TAPA + Autobridge
+ Following [Install TAPA](https://tapa.readthedocs.io/en/release/installation.html) to install TAPA(Autobridge) and Gurobi.
+ Vitis 2021.2
+ Xilinx xilinx_u280_xdma_201920_3 shell and a Xilinx U280 FPGA card.

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
    TAPAB=./bitstream/Callipepla_xilinx_u280_xdma_201920_3.xclbin ./serpenscg ../matrices/mhd3200b/mhd3200b.mtx 100

For Callipepla, 

    cd Callipepla
    TAPAB=./bitstream/Callipepla_xilinx_u280_xdma_201920_3.xclbin ./serpenscg ../matrices/mhd3200b/mhd3200b.mtx 100

