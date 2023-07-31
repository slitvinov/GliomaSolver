#!/bin/sh

: ${brain=./brain}
N=4
OMP_NUM_THREADS=$N LD_LIBRARY_PATH=VTK/lib/vtk-5.2 "$brain" \
    -adaptive 1 \
    -dumpfreq 50 \
    -Dw 0.0013 \
    -icx 0.28 \
    -icy 0.75 \
    -icz 0.35 \
    -model RD \
    -nthreads $N \
    -PatFileName ./ \
    -profiler 1 \
    -rho 0.025 \
    -Tend 300 \
    -verbose 1 \
    -vtk 1 \
