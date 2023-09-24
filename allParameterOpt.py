#!/usr/bin/python

import cmaes
import glioma_solver
import numpy as np
import os
import scipy
import struct
import sys
import nibabel as nib
import time


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
    rho = x[4]
    dw = x[3]
    period = 0
    glioma_solver.run(bpd, GM, WM, ic, dw, rho, tend, period, HG)


def fun(x):
    #print("fun run x:", x)
    sim(x)
    err = np.linalg.norm(HG - PET[::2, ::2, ::2])
    # TODO somehow 0... 
    #err = 10# diceLikelihood(HG, FLAIR[::2, ::2, ::2], T1c[::2, ::2, ::2]) + 0.00001 
    #err = np.random.rand()
    sys.stderr.write("allParameterOpt.py: %d: %.16e: %s\n" % (err, os.getpid(), str(x)))

    #print("opt.py: err:", err, "process:", os.getpid(), str(x))# %.16e: %s\n" % (err, os.getpid(), str(x)))
    return err


if __name__ == '__main__':
    print("start")
    #time
    start = time.time()
    bpd = 16
    GM = readNii("GM.nii.gz")
    WM = readNii("WM.nii.gz")
    PET = readNii("tumPET.nii.gz")
    T1c = readNii("tumT1c.nii.gz")
    FLAIR = readNii("tumFLAIR.nii.gz")
    ic0 = np.divide(scipy.ndimage.center_of_mass(PET), np.shape(PET))

    rho0 = 0.025
    dw0 = 0.03
    tend = 100

    # ranges from LMI paper with T = 100
    parameterRanges = [[0, 1], [0, 1], [0, 1], [0.0001, 0.225], [0.001, 3]] 

    HG = np.empty_like(GM, shape=(8 * bpd, 8 * bpd, 8 * bpd))
    opt = cmaes.cmaes(fun, ( *ic0, dw0, rho0), 0.05, 3, workers=0, trace=False, parameterRange=parameterRanges)
    end = time.time()	
    print("time: ", end - start)
    print("done1")
    sim(opt)
    write(HG)
    writeNii(HG)
    print("done")
