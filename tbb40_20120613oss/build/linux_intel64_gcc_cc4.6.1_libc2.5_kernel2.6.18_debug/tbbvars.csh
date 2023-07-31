#!/bin/csh
setenv TBBROOT "/cluster/home/infk/mgazzola/mattia/LIB/tbb40_20120613oss" #
setenv tbb_bin "/cluster/home/infk/mgazzola/mattia/LIB/tbb40_20120613oss/build/linux_intel64_gcc_cc4.6.1_libc2.5_kernel2.6.18_debug" #
if (! $?CPATH) then #
    setenv CPATH "${TBBROOT}/include" #
else #
    setenv CPATH "${TBBROOT}/include:$CPATH" #
endif #
if (! $?LIBRARY_PATH) then #
    setenv LIBRARY_PATH "${tbb_bin}" #
else #
    setenv LIBRARY_PATH "${tbb_bin}:$LIBRARY_PATH" #
endif #
if (! $?LD_LIBRARY_PATH) then #
    setenv LD_LIBRARY_PATH "${tbb_bin}" #
else #
    setenv LD_LIBRARY_PATH "${tbb_bin}:$LD_LIBRARY_PATH" #
endif #
 #
