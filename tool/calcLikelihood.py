#%%
import matplotlib.pyplot as plt
import numpy as np


def probabilityOfDetectedTumor(x, threshold, sigma):
    # x: true tumor concentration
    # threshold: threshold for detection
    # sigma: width of the sigmoid

    # return: probability of detection with MRI
    alpha = 0.5 + 0.5 * np.sign(x-threshold) * (1. - np.exp(- (x-threshold).astype(np.float128)**2 / sigma**2))
    return np.clip(alpha, 0.00000000000001, 0.99999999999999)
    #return np.clip( 0.5 + 0.5 * np.sign(x-threshold) * (1. - np.exp(- (x-threshold)**2 / sigma**2)), 0.00001, 0.99999 )

def logGaussianPrior(x, xPredicted, stdPredicted):
    # x (1D array): values to calculate prior for
    # xPredicted (1D array): value known with a certain unserainty as gaussien
    # stdPredicted (1D array): standard deviation of the gaussian
    # return: prior probability for x

    x = np.array(x)
    xPredicted = np.array(xPredicted)
    stdPredicted = np.array(stdPredicted)

    factor = 1 / np.sqrt(2 * np.pi * stdPredicted**2)
    return np.sum(  np.log( factor * np.exp(- (x - xPredicted)**2 / (2 * stdPredicted**2)) ))

def logLikelihood( proposedDistribution, seg, threshold, sigma):

    alphas = probabilityOfDetectedTumor(proposedDistribution, threshold, sigma)
    bernoulli = np.sum( seg.astype(float) * np.log( alphas)  + ( 1-seg) * np.log((1 - alphas))) 
    
    return bernoulli

def logPosterior(proposedDistribution, flair_seg, t1_seg, flair_threshold, t1_threshold,  sigma, addPrior = False, x = None, xMeasured = None, stdMeasured = None):
    #TODO test

    logLikelihoodFlair = logLikelihood( proposedDistribution, flair_seg, flair_threshold, sigma)
    logLikelihoodT1 = logLikelihood( proposedDistribution, t1_seg, t1_threshold, sigma)

    if addPrior:
        logPrior = logGaussianPrior(x, xMeasured, stdMeasured)
        return logLikelihoodFlair + logLikelihoodT1 +logPrior
    
    else:
        return logLikelihoodFlair + logLikelihoodT1

#dice
def dice(a, b):
    boolA, boolB = a > 0, b > 0 

    return 2 * np.sum( np.logical_and(boolA, boolB)) / (np.sum(boolA) + np.sum(boolB))

# in addition it might make sense to use the minumum dice for several thresholds instead of infering the thresholds
def diceLogLikelihood(proposedDistribution, flair_seg, t1_seg):
    # TODO not working
    th_flair = 0.3
    th_core = 0.75
    diceFlair = dice(proposedDistribution > th_flair, flair_seg) 
    diceT1 = dice(proposedDistribution > th_core, t1_seg)

    # add this to start convergence and prevent 0 dice
    if diceFlair +diceT1 < 0.001:
        additionalStartingDice = np.log(0.1 * dice(proposedDistribution > 0.01, flair_seg) + 0.0000001)
    else:
        additionalStartingDice = 0

    print("additionalStartingDice:", additionalStartingDice)

    return np.log(diceFlair+ 0.0000001) + np.log(diceT1+ 0.0000001) + additionalStartingDice


# %% check dice and probabiltiy of detectioy  
if __name__ == '__main__':
    def readNii(path):
        return nib.load(path).get_fdata().astype(np.float32)
    import nibabel as nib

    #path = "results/2023_09_26-01_35_51_gen_1000/result.nii.gz"
    tumor = nib.load("/home/jonas/workspace/programs/GliomaSolver/results/2023_09_26-01_35_51_gen_1000/result.nii.gz").get_fdata().astype(np.float32)
    GM = readNii("../GM.nii.gz")
    WM = readNii("../WM.nii.gz")
    T1c = readNii("../tumT1c.nii.gz")
    FLAIR = readNii("../tumFLAIR.nii.gz")

    diceFlair = dice(tumor > 0.3, FLAIR[::2, ::2, ::2])
    diceT1 = dice(tumor > 0.75, T1c[::2, ::2, ::2])

    tumor = 0 * tumor+1
    #print(logLikelihood(tumor, T1c[::2, ::2, ::2], 0.3, 0.05))
    err = -logPosterior(tumor, FLAIR[::2, ::2, ::2], T1c[::2, ::2, ::2], 0.3, 0.75, 0.05, addPrior = False)
    errFlair = -logLikelihood(tumor, FLAIR[::2, ::2, ::2], 0.3, 0.05)
    errT1 = -logLikelihood(tumor, T1c[::2, ::2, ::2], 0.75, 0.05)
    print('err=',"%.2e" % errFlair, "%.2e" % errT1)
    print('errTot=',"%.2e" % err)
    errZero = -logPosterior(np.zeros_like(tumor), FLAIR[::2, ::2, ::2], T1c[::2, ::2, ::2], 0.3, 0.75, 0.05, addPrior = False)
    print('errZero=',"%.2e" % errZero)
    errFlairFlair = -logLikelihood(FLAIR[::2, ::2, ::2], FLAIR[::2, ::2, ::2], 0.3, 0.05)
    errT1Flair = -logLikelihood(FLAIR[::2, ::2, ::2], T1c[::2, ::2, ::2], 0.75, 0.05)
    print('errFlair=',"%.2e" % errFlairFlair, "%.2e" % errT1Flair)
    #%%
    #readNii = lambda path: nib.load(path).get_fdata().astype(np.float32)
    #posterior(np.linspace(0,1,100), np.ones(100), np.ones(100), 0.25, 0.75, 0.05)
    np.log(probabilityOfDetectedTumor(np.array([0]), 0.3, 0.05))
    # %%
# %% 
if __name__ == '__main__':
    # janas paper example
    x = np.linspace(0,1,100)
    titles = ["FLAIR", "T1c"]
    thresholds = [0.25, 0.75]
    for i in range(2):
        y = probabilityOfDetectedTumor(x,thresholds[i],0.08)
        plt.title("Jana - probability of detected tumor " )
        plt.plot(x,y, label = titles[i])
        plt.ylabel("Probability of Detection") 
        plt.xlabel("True Tumor Concentration")

    plt.legend()
# %%
