#!/bin/bash

cmake -S . -B build
cmake --build build --target rmsd
#cmake --build build
