#!/usr/bin/python

import cmaes
import glioma_solver
import numpy as np
import os
import scipy
import sys
import nibabel as nib


def readNii(path):
    return nib.load(path).get_fdata().astype(np.float32)


def writeNii(array):
    path = "%dx%dx%dle.nii.gz" % np.shape(array)
    nibImg = nib.Nifti1Image(array, np.eye(4))
    nib.save(nibImg, path)


def clip(x):
    if x < 0: return x
    if x > 1: return 1
    return x


def sim(x):
    ix, iy, iz, rho, dw = x
    ix = clip(ix)
    iy = clip(iy)
    iz = clip(iz)
    dw = abs(dw)
    rho = abs(rho)
    glioma_solver.run(bpd, GM, WM, (ix, iy, iz), dw, rho, tend, period, HG)


def fun(x):
    sim(x)
    err = np.linalg.norm(HG - PET[::2, ::2, ::2])
    sys.stderr.write("opt1.py: %d: %.16e: %s\n" % (os.getpid(), err, str(x)))
    return err


bpd = 16
period = 0
GM = readNii("GM.nii.gz")
WM = readNii("WM.nii.gz")
PET = readNii("tumPET.nii.gz")
ic0 = np.divide(scipy.ndimage.center_of_mass(PET), np.shape(PET))
rho0 = 0.025
dw0 = 0.0013
tend = 300
HG = np.empty_like(GM, shape=(8 * bpd, 8 * bpd, 8 * bpd))
opt = cmaes.cmaes(fun, (*ic0, rho0, dw0), sigma=0.05, g_max=100, workers=3)
period = 10
sim(opt)
writeNii(HG)
