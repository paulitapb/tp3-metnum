#!/bin/bash

rm -rf *.o
cmake -DPYTHON_EXECUTABLE="$(which python3)" -DCMAKE_BUILD_TYPE=Release ..
make
echo "FINISHED COMPILATION"
echo "Run ./tp3 to run tests"