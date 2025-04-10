#!/bin/bash

make runclean
make clean
rm -rf CMakeFiles CMakeCache.txt debug results
mkdir debug results
cmake .
make clean
make release
