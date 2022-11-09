#!/bin/bash

cd build
rm -rf *
cmake -DPYTHON_EXECUTABLE="$(which python3)" -DCMAKE_BUILD_TYPE=Release ..
make
mv ?*.so ../