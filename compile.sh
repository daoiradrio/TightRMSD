#!/bin/bash

cmake -S . -B build
cmake --build build --target tight-rmsd
#cmake --build build
