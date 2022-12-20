MAT_PATH='/share/slh/CG_DATA/matrices_serpenscg'
XCLBIN_PATH='../SerpensCG/bitstream/SerpensCG_xilinx_u280_xdma_201920_3.xclbin'
g++ -o serpenscg -Wno-write-strings -Wunused-result -O2 ../SerpensCG/src/serpenscg.cpp ./serpenscg-host.cpp -ltapa -lfrt -lglog -lgflags -lOpenCL

TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/ex9.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/bcsstk15.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/bodyy4.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/ted_B.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/ted_B_unscaled.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/bcsstk24.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/nasa2910.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/s3rmt3m3.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/bcsstk28.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/s2rmq4m1.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/cbuckle.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/olafu.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/gyro_k.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/bcsstk36.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/msc10848.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/raefsky4.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/nd3k.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/nd6k.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/2cubes_sphere.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/cfd2.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/Dubcova3.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/ship_003.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/offshore.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/shipsec5.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/ecology2.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/tmt_sym.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/boneS01.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/hood.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/bmwcra_1.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/af_shell3.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/Fault_639.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/Emilia_923.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/Geo_1438.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/Serena.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/audikw_1.bin 20000
TAPAB=${XCLBIN_PATH} ./serpenscg ${MAT_PATH}/Flan_1565.bin 20000
