# CMSC858D HW1

## Overview
This HW consists of implementations for rank support on a bit vector, select support on the same vector, and a sparse array. The implementations build on each other and are found in their respective directories. In addition, hw1.cpp runs/tests all of them.

## Directions to Run
- Make sure SDSL is installed (the include/ and lib/ directories are accessible)
- Run below commands to compile and then run the binary (can use g++ instead of clang)
  - Replace "./include" and "./lib" with the SDSL include and lib paths on your system

  - `clang++ -std=c++20 -O3 -DNDEBUG -I ./include -L ./lib hw1.cpp -o hw1 -lsdsl -ldivsufsort -ldivsufsort64`
  - `./hw1`

- Code will output whether each class works correctly and writes to the testing files (used for plots)
