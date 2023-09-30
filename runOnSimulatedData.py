
#%%
from  recPred.dataset_syntheticBrats import Dataset_SyntheticBrats
from allParameterOpt import readNii, writeNii, CmaesSolver
import numpy as np
import os 
import time
#%%
def extendTo256(a):
    return np.repeat(np.repeat(np.repeat(a, 2, axis=0), 2, axis=1), 2, axis=2)

if __name__ == '__main__':
    print("start")


    directoryPath = "/mnt/8tb_slot8/jonas/datasets/lmiSynthetic/"

    testDataset = Dataset_SyntheticBrats(directoryPath, "/mnt/8tb_slot8/jonas/workingDirDatasets/lmiSynthetic", deleteAndReloadWorkingDir=True)
    patientID = 0

    GM128 = testDataset.loadPatientImageEnumarated(patientID, "GM", '0', "dat128JanasSolver").astype(np.float32)
    GM = GM128 #GM = np.asfortranarray(extendTo256(GM128))

    WM128 = testDataset.loadPatientImageEnumarated(patientID, "WM", '0', "dat128JanasSolver").astype(np.float32)
    WM = WM128 #np.asfortranarray(extendTo256(WM128))

    T1c = testDataset.loadPatientImageEnumarated(patientID, "seg-t1c", '0', "dat128JanasSolver")[:-1,:-1,:-1]
    FLAIR = testDataset.loadPatientImageEnumarated(patientID, "seg-flair", '0', "dat128JanasSolver")[:-1,:-1,:-1]



    settings = {}
    # ranges from LMI paper with T = 100
    parameterRanges = [[0, 1], [0, 1], [0, 1], [0.0001, 0.225], [0.001, 3]] 
    settings["parameterRanges"] = parameterRanges

    settings["bpd"] = 16
    settings["rho0"] = 0.1
    settings["dw0"] = 0.2
    settings["workers"] = 2
    settings["sigma0"] = 0.02
    settings["generations"] = 500
    settings["lossfunction"] = "bernoulli"#"dice"#

    debugMode = False
    if debugMode == True:
        settings["generations"] = 1
        settings["workers"] = 0

    settings["addPrior"] = False
    
    #TODO delet None in getParams
    params = testDataset.getParams(patientID, None)

    # the factor to increase the range for the ensemble
    settings["factorSTD"] = 3
    settings["stdMeasured"], settings["xMeasured"] =[], []
    for key in ["x", "y", "z", "D-cm-d", "rho1-d"]:
        print(key, params[key]['ensemble'])

        settings["xMeasured"].append(np.mean(params[key]['ensemble']))
        settings["stdMeasured"].append(np.std(params[key]['ensemble']) * settings["factorSTD"])

    if settings["addPrior"]:
        settings["rho0"] = settings["xMeasured"][4]
        settings["dw0"] = settings["xMeasured"][3]

    solver = CmaesSolver(settings, WM, GM, FLAIR, T1c)
    resultTumor, resultDict = solver.run()

    # save results
    datetime = time.strftime("%Y_%m_%d-%H_%M_%S")
    path = "./resultsSynData/"+ 'patient_' + str(patientID) + '_dtime' + datetime +"_gen_"+ str(settings["generations"]) + "_loss_" + str(settings["lossfunction"]) + '_prior_' + str(settings["addPrior"]) + "_dw0_" + str(settings["dw0"]) + "_rho0_" + str(settings["rho0"])+"/"
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

# %%
