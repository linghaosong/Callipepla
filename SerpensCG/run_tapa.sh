tapac \
  --work-dir run \
  --top SerpensCG \
  --part-num xcu280-fsvh2892-2L-e \
  --platform xilinx_u280_xdma_201920_3 \
  --clock-period 4.00 \
  -o SerpensCG.xo \
  --constraint SerpensCG_floorplan.tcl \
  --connectivity ../link_config.ini \
  --read-only-args edge_list_ptr \
  --read-only-args edge_list_ch* \
  --read-only-args vec_digA \
  --write-only-args vec_res \
  --enable-synth-util \
  --max-parallel-synth-jobs 16 \
  --enable-hbm-binding-adjustment \
  --run-floorplan-dse \
  --min-area-limit 0.58 \
  --min-slr-width-limit 5000 \
  ../src/serpenscg.cpp \
   2>&1 | tee tapa.log
