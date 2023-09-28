
#%%
from  recPred.dataset_syntheticBrats import Dataset_SyntheticBrats
from allParameterOpt import readNii, writeNii, CmaesSolver
import numpy as np
import os 
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
    settings["generations"] =50
    settings["lossfunction"] = "bernoulli"#"dice"#
    print('Lossfunction:', settings["lossfunction"])

    solver = CmaesSolver(settings, WM, GM, FLAIR, T1c)
    resultTumor, resultDict = solver.run()

    # save results
    datetime = time.strftime("%Y_%m_%d-%H_%M_%S")
    path = "./resultsSynData/"+ datetime +"_gen_"+ str(settings["generations"]) + "_loss_" + str(settings["lossfunction"]) + "/"
    os.makedirs(path, exist_ok=True)
    np.save(path + "settings.npy", settings)
    np.save(path + "results.npy", resultDict)
    writeNii(resultTumor, path = path+"result.nii.gz")
    
    print("diceT1_67",  resultDict["diceT1_67"])
    print("diceFLAIR_25",  resultDict["diceFLAIR_25"])
    print("likelihoodFlair_25",  resultDict["likelihoodFlair_25"])
    print("likelihoodT1_75",  resultDict["likelihoodT1_75"])
    print("final_loss",  resultDict["final_loss"])
    print("opt_params",  resultDict["opt_params"])


#%% test

if __name__ == '__main__':

    directoryPath = "/mnt/8tb_slot8/jonas/datasets/lmiSynthetic/"

    testDataset = Dataset_SyntheticBrats(directoryPath, "/mnt/8tb_slot8/jonas/workingDirDatasets/lmiSynthetic", deleteAndReloadWorkingDir=True)
    patientID = 2
    gt = testDataset.loadPatientImageEnumarated(patientID, "groundTruth", '0', "dat128JanasSolver")
    wm = testDataset.loadPatientImageEnumarated(patientID, "WM", '0', "dat128JanasSolver")

    csf = testDataset.loadPatientImageEnumarated(patientID, "CSF", '0', "dat128JanasSolver")
    gm = testDataset.loadPatientImageEnumarated(patientID, "GM", '0', "dat128JanasSolver")

    seg_flair = testDataset.loadPatientImageEnumarated(patientID, "seg-flair", '0', "dat128JanasSolver")

    seg_t1c = testDataset.loadPatientImageEnumarated(patientID, "seg-t1c", '0', "dat128JanasSolver")
    
    import matplotlib.pyplot as plt
    from scipy import ndimage
    centerofmass = np.array(ndimage.measurements.center_of_mass(gt))
    slice = int(centerofmass[2]) 
    plt.imshow(wm[:,:,slice], cmap='gray')
    
    plt.imshow(gt[:,:,slice], alpha=gt[:,:,slice] , cmap='hot')
    plt.imshow(seg_t1c[:,:,slice], alpha=0.2 , cmap='gray')
    plt.imshow(seg_flair[:,:,slice], alpha=0.2  , cmap='gray')


# %%

# %%
import pkgutil

package = recPred
for importer, modname, ispkg in pkgutil.iter_modules(package.__path__):
    print(modname)
# %%
