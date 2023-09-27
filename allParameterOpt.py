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


class CmaesSolver():
    def __init__(self,  settings, wm, gm, flair, t1c):
        lossfunction = settings["lossfunction"]
        if lossfunction == "dice":
            self.lossfunction = calcLikelihood.diceLogLikelihood
        elif lossfunction == "bernoulli":
            self.lossfunction = calcLikelihood.logPosterior

        self.tend = 100
        self.flair_th = 0.25
        self.t1c_th = 0.675

        self.T1c = t1c
        self.FLAIR = flair
        self.GM = gm
        self.WM = wm

        self.settings = settings

        self.bpd = settings["bpd"]
        
        self.rho0 = settings["rho0"]
        self.dw0 = settings["dw0"]
        self.workers = settings["workers"]
        self.sigma0 = settings["sigma0"]
        self.generations = settings["generations"]
        self.parameterRanges = settings["parameterRanges"]
        self.HG = np.empty_like(self.GM, shape=(8 * self.bpd, 8 * self.bpd, 8 * self.bpd))

    def sim(self, x):
        ic = x[:3]
        rho = x[4]
        dw = x[3]
        period = 0
        glioma_solver.run(self.bpd, self.GM, self.WM, ic, dw, rho, self.tend, period, self.HG)


    def fun(self, x):
        self.sim(x)
        err =  - self.lossfunction(self.HG, self.FLAIR, self.T1c, self.flair_th, self.t1c_th) 
        
        # TODO compare original with dice likelihood on synthetic data
        sys.stderr.write("allParameterOpt.py: %d: %.16e: %s\n" % (err, os.getpid(), str(x)))
        print(err, " --  run err:")

        return err

    def run(self):
        start = time.time()

        # set initial conditions
        ic0 = np.divide(scipy.ndimage.center_of_mass(self.FLAIR), np.shape(self.FLAIR))
        
        opt = cmaes.cmaes(self.fun, ( *ic0, self.dw0, self.rho0), self.sigma0, self.generations, workers=self.workers, trace=False, parameterRange= self.parameterRanges)
        self.sim(opt)
        end = time.time()

        resultDict = {}

        resultDict["diceT1_75"] = calcLikelihood.dice(self.HG > 0.675, self.T1c)
        resultDict["diceFLAIR_30"] = calcLikelihood.dice(self.HG > 0.25, self.FLAIR)
        resultDict["likelihoodFlair_30"] = calcLikelihood.logLikelihood(self.HG, self.FLAIR, 0.25, 0.05)
        resultDict["likelihoodT1_75"] = calcLikelihood.logLikelihood(self.HG, self.T1c, 0.675, 0.05)
        resultDict['final_loss'] = self.fun(opt)
        
        resultDict["opt_params"] = opt
        resultDict["time_min"] = (end - start) / 60
        
        return self.HG, resultDict

        
#%%
if __name__ == '__main__':
    print("start")

    GM = readNii("GM.nii.gz")#[::2, ::2, ::2]
    WM = readNii("WM.nii.gz")#[::2, ::2, ::2]
    T1c = readNii("tumT1c.nii.gz")[::2, ::2, ::2]
    FLAIR = readNii("tumFLAIR.nii.gz")[::2, ::2, ::2]

    settings = {}
    # ranges from LMI paper with T = 100
    parameterRanges = [[0, 1], [0, 1], [0, 1], [0.0001, 0.225], [0.001, 3]] 
    settings["parameterRanges"] = parameterRanges

    settings["bpd"] = 16
    settings["rho0"] = 0.025
    settings["dw0"] = 0.2
    settings["workers"] = 8
    settings["sigma0"] = 0.05
    settings["generations"] =1000
    settings["lossfunction"] = "bernoulli"#"dice"#
    print('Lossfunction:', settings["lossfunction"])

    solver = CmaesSolver(settings, WM, GM, FLAIR, T1c)
    resultTumor, resultDict = solver.run()

    # save results
    datetime = time.strftime("%Y_%m_%d-%H_%M_%S")
    path = "./results/"+ datetime +"_gen_"+ str(settings["generations"]) + "_loss_" + str(settings["lossfunction"]) + "/"
    os.makedirs(path, exist_ok=True)
    np.save(path + "settings.npy", settings)
    np.save(path + "results.npy", resultDict)
    writeNii(resultTumor, path = path+"result.nii.gz")
    
    print("diceT1_75",  resultDict["diceT1_75"])
    print("diceFLAIR_30",  resultDict["diceFLAIR_30"])
    print("likelihoodFlair_30",  resultDict["likelihoodFlair_30"])
    print("likelihoodT1_75",  resultDict["likelihoodT1_75"])
    print("final_loss",  resultDict["final_loss"])
    print("opt_params",  resultDict["opt_params"])

#%%
if False: #__name__ == '__main__':
    dict = np.load("/home/jonas/workspace/programs/GliomaSolver/results/2023_09_26-17_11_36_gen_10/results.npy", allow_pickle=True).item()
    params = dict['opt_params']
    #solver = CmaesSolver(settings, WM, GM, FLAIR, T1c, likelihoodFunction = "dice")
    #arrray = solver.sim(params)

    writeNii(WM[::2,::2,::2], path = "./results/WM.nii.gz")
    writeNii(GM[::2,::2,::2], path = "./results/GM.nii.gz")
    writeNii(T1c, path = "./results/tumT1c.nii.gz")
    writeNii(FLAIR, path = "./results/tumFLAIR.nii.gz")
    # %%

# %%
