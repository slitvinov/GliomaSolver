#!/usr/bin/python

import sys
import struct
import glioma_solver
import numpy as np


def unpack(string, file):
    buffer = file.read(struct.calcsize(string))
    return struct.unpack(string, buffer)


def read(path):
    with open(path, "rb") as inp:
        magic, = unpack('i', inp)
        if magic != 1234:
            sys.stderr.write("fun.py error: '%s' not a dat file\n" % path)
            sys.exit(1)
        dim, = unpack('i', inp)
        if dim != 3:
            sys.stderr.write("fun.py error: wrong dimension '%d'\n" % dim)
            sys.exit(1)
        nx, ny, nz, type_id = unpack('iiii', inp)
        seek = inp.tell()
        type_name = {0: "<d", 1: "<f", 2: "<i"}[type_id]
    with open(path, "rb+") as inp:
        mm = inp.read()
        dtype = np.dtype(type_name)
        return np.ndarray((nx, ny, nz), dtype, mm, seek, order='F')


def write(a):
    path = "%dx%dx%dle.raw" % np.shape(a)
    a.tofile(path)
    print(path)


bpd = 32
GM = read("GM.dat")
WM = read("WM.dat")
HG = np.empty_like(GM, shape=(8 * bpd, 8 * bpd, 8 * bpd))
ic = 0.64925073, 0.59093041, 0.37127833
dw = 0.0013
rho = 0.025
tend = 300
glioma_solver.run(bpd, GM, WM, ic, dw, rho, tend, HG)

PET = read("tumPET.dat")
T1c = read("tumT1c.dat")
FLAIR = read("tumFLAIR.dat")

PETsigma2 = 0.000361
PETscale = 0.85
T1uc = 0.7
T2uc = 0.25
slope = 2
ans = glioma_solver.likelihood(HG, PET, T1c, FLAIR, PETsigma2, PETscale, slope,
                               T1uc, T2uc)
print("%.16e" % ans)
write(HG)
