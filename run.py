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
            sys.stderr.write("run.py: error: '%s' not a dat file\n" % path)
            sys.exit(1)
        dim, = unpack('i', inp)
        if dim != 3:
            sys.stderr.write("run.py: error: wrong dimension '%d'\n" % dim)
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

GM = read("GM.dat")
WM = read("WM.dat")
HG = np.empty_like(GM, shape=(128, 128, 128))
ic = 0.6497946102507519, 0.5908331665234543, 0.3715947899171972
dw = 0.0013
rho = 0.025
tend = 300

glioma_solver.run(GM, WM, ic, dw, rho, tend, HG)


write(HG)

