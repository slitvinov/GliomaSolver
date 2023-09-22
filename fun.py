#!/usr/bin/python
#%%
import sys
import struct
import glioma_solver
import numpy as np
import nibabel as nib


def unpack(string, file):
    buffer = file.read(struct.calcsize(string))
    return struct.unpack(string, buffer)


def readDat(path):
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
    
def readNii(path):
    return nib.load(path).get_fdata().astype(np.float32)

def write(a):
    path = "%dx%dx%dle.raw" % np.shape(a)
    a.tofile(path)
    print(path)

def writeNii(array):
    path = "%dx%dx%dle.nii.gz" % np.shape(array)
    nibImg = nib.Nifti1Image(array, np.eye(4))
    nib.save(nibImg, path)

bpd = 32
GM = readNii('GM.nii.gz')
WM = readNii('WM.nii.gz')
HG = np.empty_like(GM, shape=(8 * bpd, 8 * bpd, 8 * bpd))
ic = 0.64925073, 0.59093041, 0.37127833
dw = 0.0013
rho = 0.025
tend = 300
period = 10
glioma_solver.run(bpd, GM, WM, ic, dw, rho, tend, period, HG)

PET = readNii("tumPET.nii.gz")
T1c = readNii("tumT1c.nii.gz")
FLAIR = readNii("tumFLAIR.nii.gz")

PETsigma2 = 0.000361
PETscale = 0.85
T1uc = 0.7
T2uc = 0.25
slope = 2
ans = glioma_solver.likelihood(HG, PET, T1c, FLAIR, PETsigma2, PETscale, slope,
                               T1uc, T2uc)
print("%.16e" % ans)
write(HG)
writeNii(HG)
np.save("HG.npy", HG)

# %%
