#%%
import numpy as np  
import matplotlib.pyplot as plt
import scipy.ndimage
import opt
import nibabel as nib

hg = np.load('HG.npy')


#center of mass
com = scipy.ndimage.center_of_mass(hg), np.shape(hg)
plt.imshow(hg[166,:,:])
# %%
testDatRead = opt.read('/home/jonas/workspace/programs/GliomaSolver/testFilesJonas/tgm002_dat/WM.dat') 

testNii = nib.load('/home/jonas/workspace/programs/GliomaSolver/testFilesJonas/tgm002/WM.nii').get_fdata()

# %%
slice = 80
plt.imshow(testDatRead[:,:,slice])
plt.show()
plt.imshow(testNii[:,:,slice])
# %%
plt.imshow((testNii - np.roll(testDatRead, 0, axis=0))[:,:,slice])
plt.colorbar()
# %%
plt.imshow((testNii - np.flip(testDatRead, axis=1))[:,:,slice])

# %% just read dat files and save as nifti
wmDat = opt.read('WM.dat') 
gmDat = opt.read('GM.dat')

PETDat = opt.read("tumPET.dat")
T1cDat = opt.read("tumT1c.dat")
FLAIRDat = opt.read("tumFLAIR.dat")


# save nifti 
wmNii = nib.Nifti1Image(wmDat, np.eye(4))
gmNii = nib.Nifti1Image(gmDat, np.eye(4))
PETNii = nib.Nifti1Image(PETDat, np.eye(4))
T1cNii = nib.Nifti1Image(T1cDat, np.eye(4))
FLAIRNii = nib.Nifti1Image(FLAIRDat, np.eye(4))


nib.save(wmNii, 'WM.nii.gz')
nib.save(gmNii, 'GM.nii.gz')
nib.save(PETNii, 'tumPET.nii.gz')
nib.save(T1cNii, 'tumT1c.nii.gz')
nib.save(FLAIRNii, 'tumFLAIR.nii.gz')

# %%
print(np.max(nib.load('WM.nii.gz').get_fdata() - wmDat))
print(np.max(nib.load('GM.nii.gz').get_fdata() - gmDat))
print(np.max(nib.load('tumPET.nii.gz').get_fdata() - PETDat))
print(np.max(nib.load('tumT1c.nii.gz').get_fdata() - T1cDat))
print(np.max(nib.load('tumFLAIR.nii.gz').get_fdata() - FLAIRDat))

# %%
