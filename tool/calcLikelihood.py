#%%
import matplotlib.pyplot as plt
import numpy as np


def probabilityOfDetectedTumor(x, threshold, sigma):
    # x: true tumor concentration
    # threshold: threshold for detection
    # sigma: width of the sigmoid

    # return: probability of detection with MRI
    return 0.5 + 0.5 * np.sign(x-threshold) * (1. - np.exp(- (x-threshold)**2 / sigma**2))

def gaussianPrior(x, xPredicted, stdPredicted):
    # x (1D array): values to calculate prior for
    # xPredicted (1D array): value known with a certain unserainty as gaussien
    # stdPredicted (1D array): standard deviation of the gaussian
    # return: prior probability for x

    # TODO implement
    return 1

def likelihood(x, xMeasured, stdMeasured):
    #TODO   
    return 1

#dice
def dice(a, b):
    return 2 * np.sum(a * b) / (np.sum(a) + np.sum(b))

# in addition it might make sense to use the minumum dice for several thresholds instead of infering the thresholds
def diceLikelihood(proposedDistribution, flair_seg, t1_seg):
    # TODO not working
    th_flair = 0.3
    th_core = 0.75
    diceFlair = dice(proposedDistribution > th_flair, flair_seg)
    diceT1 = dice(proposedDistribution > th_core, t1_seg)
    return diceFlair * diceT1

# %% 
if __name__ == '__main__':
    # janas paper example
    x = np.linspace(0,1,100)
    titles = ["FLAIR", "T1c"]
    thresholds = [0.25, 0.75]
    for i in range(2):
        y = probabilityOfDetectedTumor(x,thresholds[i],0.05)
        plt.title("Jana - probability of detected tumor " )
        plt.plot(x,y, label = titles[i])
        plt.ylabel("Probability of Detection") 
        plt.xlabel("True Tumor Concentration")

    plt.legend()
# %%
