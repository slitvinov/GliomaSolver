#!/bin/sh

: ${brain=./brain}

"$brain" \
    -adaptive 1 \
    -dumpfreq 50 \
    -model RD \
    -nthreads $N \
    -PatFileName ./ \
    -profiler 1 \
    -UQ 1 \
    -verbose 1 \
    -vtk 1 \
