#!/usr/bin/python
#%%
import cmaes
import glioma_solver
import numpy as np
import os
import scipy
import struct
import sys
import nibabel as nib
import time
from tool import calcLikelihood


def readNii(path):
    return nib.load(path).get_fdata().astype(np.float32)

def write(a):
    path = "%dx%dx%dle.raw" % np.shape(a)
    with open(path, "wb") as file:
        file.write(a.tobytes('F'))
    sys.stderr.write("opt.py: write: %s\n" % path)

def writeNii(array, path = ""):
    if path == "":
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
    #err = np.linalg.norm(HG - PET[::2, ::2, ::2])
    err =  - calcLikelihood.diceLogLikelihood(HG, FLAIR[::2, ::2, ::2], T1c[::2, ::2, ::2]) 
    
    # TODO compare original with dice likelihood on synthetic data
    sys.stderr.write("allParameterOpt.py: %d: %.16e: %s\n" % (err, os.getpid(), str(x)))
    print(err, " --  run err:")

    return err

#%%
if __name__ == '__main__':
    settings = {}
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
    dw0 = 0.2
    tend = 100
    workers = 8
    sigma0 = 0.05
    generations = 1000

    # ranges from LMI paper with T = 100
    parameterRanges = [[0, 1], [0, 1], [0, 1], [0.0001, 0.225], [0.001, 3]] 
    
    HG = np.empty_like(GM, shape=(8 * bpd, 8 * bpd, 8 * bpd))
    opt = cmaes.cmaes(fun, ( *ic0, dw0, rho0), sigma0, generations, workers=workers, trace=False, parameterRange=parameterRanges)
    end = time.time()	
    print("done time: ", (end - start) / 60, "min") 
    print("parameters:", opt)
    sim(opt)
    write(HG)

    # save results
    datetime = time.strftime("%Y_%m_%d-%H_%M_%S")
    path = "./results/"+ datetime +"_gen_"+ str(generations) + "/"

    os.makedirs(path, exist_ok=True)
    
    
    print("dice T1:", calcLikelihood.dice(HG > 0.75, T1c[::2, ::2, ::2]))
    print("dice FLAIR:", calcLikelihood.dice(HG > 0.3, FLAIR[::2, ::2, ::2]))
    settings["diceT1"] = calcLikelihood.dice(HG > 0.75, T1c[::2, ::2, ::2])
    settings["diceFLAIR"] = calcLikelihood.dice(HG > 0.3, FLAIR[::2, ::2, ::2])
    settings["opt_params"] = opt
    settings["time_min"] = (end - start) / 60
    settings["parameterRanges"] = parameterRanges
    settings["bpd"] = bpd
    settings["rho0"] = rho0
    settings["dw0"] = dw0
    settings["tend0"] = tend
    settings["workers"] = workers
    settings["sigma0"] = sigma0
    settings["generations"] = generations
    settings["initialParams"] = ( *ic0, dw0, rho0)

    np.save(path + "settings.npy", settings)
    print("saved settings: ", settings )

    writeNii(HG, path = path+"result.nii.gz")

    print("done")
