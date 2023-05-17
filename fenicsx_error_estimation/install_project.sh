#!/bin/bash

# Build c++ backend
cmake -G Ninja -B build-dir -DCMAKE_BUILD_TYPE=Release cpp/
ninja -C build-dir install

# Link c++ and python
pip3 install python/ -v --upgrade