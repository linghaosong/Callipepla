g++ -o serpenscg -Wno-write-strings -Wunused-result -O2 ../SerpensCG/src/serpenscg.cpp ./serpenscg-host.cpp -ltapa -lfrt -lglog -lgflags -lOpenCL
g++ -o callipepla -Wno-write-strings -Wunused-result -O2 ../Callipepla/src/callipepla.cpp ./callipepla-host.cpp -ltapa -lfrt -lglog -lgflags -lOpenCL
