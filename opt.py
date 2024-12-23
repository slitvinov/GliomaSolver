#!/usr/bin/python

import glioma_solver
import numpy as np
import os
import scipy
import sys
import nibabel as nib

def readNii(path):
    return nib.load(path).get_fdata().astype(np.float32)

def write(a):
    path = "%dx%dx%dle.raw" % np.shape(a)
    with open(path, "wb") as file:
        file.write(a.tobytes('F'))
    sys.stderr.write("opt.py: write: %s\n" % path)

def writeNii(array):
    path = "%dx%dx%dle.nii.gz" % np.shape(array)
    nibImg = nib.Nifti1Image(array, np.eye(4))
    nib.save(nibImg, path)

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
    GM = readNii("GM.nii.gz")
    WM = readNii("WM.nii.gz")
    PET = readNii("tumPET.nii.gz")
    ic0 = np.divide(scipy.ndimage.center_of_mass(PET), np.shape(PET))
    rho0 = 0.025
    dw = 0.0013
    tend = 300
    HG = np.empty_like(GM, shape=(8 * bpd, 8 * bpd, 8 * bpd))
    opt = scipy.optimize.differential_evolution(fun, ((0, 1), (0, 1), (0, 1),
                                                      (0.01, 0.04)),
                                                updating='deferred',
                                                x0=(*ic0, rho0),
                                                polish=False,
                                                disp=True,
                                                workers=-1,
                                                maxiter=4)
    sim(opt.x)
    write(HG)
    writeNii(HG)
