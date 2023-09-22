#!/usr/bin/python

import cmaes
import glioma_solver
import numpy as np
import os
import scipy
import struct
import sys


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
    with open(path, "wb") as file:
        file.write(a.tobytes('F'))
    sys.stderr.write("opt.py: write: %s\n" % path)


def sim(x):
    ic = x[:3]
    rho = x[3]
    period = 0
    glioma_solver.run(bpd, GM, WM, ic, dw, rho, tend, period, HG)


def fun(x):
    sim(x)
    err = np.linalg.norm(HG - PET[::2, ::2, ::2])
    sys.stderr.write("opt.py: %d: %.16e: %s\n" % (err, os.getpid(), str(x)))
    return err


if __name__ == '__main__':
    bpd = 16
    GM = read("GM.dat")
    WM = read("WM.dat")
    PET = read("tumPET.dat")
    ic0 = np.divide(scipy.ndimage.center_of_mass(PET), np.shape(PET))
    rho0 = 0.025
    dw = 0.0013
    tend = 300
    HG = np.empty_like(GM, shape=(8 * bpd, 8 * bpd, 8 * bpd))
    opt = cmaes.cmaes(fun, (*ic0, rho0), 0.05, 100, workers=2)
    sim(opt)
    write(HG)