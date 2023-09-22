#!/usr/bin/python

import sys
import matplotlib.pylab as plt
import array
import numpy as np

def read(path):
    with open(path, "rb") as inp:
        a = array.array('i')
        a.fromfile(inp, 6)
        magic, dim, nx, ny, nz, type_id = a
        if magic != 1234:
            sys.stderr.write("mip.py: error: '%s' not a dat file\n" % path)
            sys.exit(1)
        if dim != 3:
            sys.stderr.write("mip.py: error: wrong dimension '%d'\n" % dim)
            sys.exit(1)
        seek = inp.tell()
        type_name = {0: "<d", 1: "<f", 2: "<i"}[type_id]
    with open(path, "rb+") as inp:
        mm = inp.read()
        dtype = np.dtype(type_name)
        return np.ndarray((nx, ny, nz), dtype, mm, seek, order='F')

a = read(sys.argv[1])
a = np.max(a, axis=2)
plt.imshow(a.T)
plt.savefig("a.png")
