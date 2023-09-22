#!/usr/bin/python

import sys
import array
import numpy as np
import skimage.measure


def read(path):
    with open(path, "rb") as inp:
        a = array.array('i')
        a.fromfile(inp, 6)
        magic, dim, nx, ny, nz, type_id = a
        if magic != 1234:
            sys.stderr.write("surface.py: error: '%s' not a dat file\n" % path)
            sys.exit(1)
        if dim != 3:
            sys.stderr.write("surface.py: error: wrong dimension '%d'\n" % dim)
            sys.exit(1)
        seek = inp.tell()
        type_name = {0: "<d", 1: "<f", 2: "<i"}[type_id]
    with open(path, "rb+") as inp:
        mm = inp.read()
        dtype = np.dtype(type_name)
        return np.ndarray((nx, ny, nz), dtype, mm, seek, order='F')


a = read(sys.argv[1])
sx, sy, sz = 1, 1, 1
nx, ny, nz = np.shape(a)
dx = 1.0 / nx / sx
dy = 1.0 / ny / sy
dz = 1.0 / nz / sz
verts, faces, normals, values = skimage.measure.marching_cubes(a,
                                                               0.5,
                                                               spacing=(sx, sy,
                                                                        sz))
