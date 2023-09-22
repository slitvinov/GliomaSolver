#!/usr/bin/python
#%%
import sys
import glioma_solver
import numpy as np
import nibabel as nib

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
