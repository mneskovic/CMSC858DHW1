# CMSC858D HW1

## Directions to Run
- Make sure SDSL is installed (the include/ and lib/ directories are accessible)
- Run below commands to compile and then run the binary (can use g++ instead of clang)
  - Replace "./include" and "./lib" to the SDSL include and lib paths on your system

  - `clang++ -std=c++20 -O3 -DNDEBUG -I ./include -L ./lib hw1.cpp -o hw1 -lsdsl -ldivsufsort -ldivsufsort64`
  - `./hw1`
