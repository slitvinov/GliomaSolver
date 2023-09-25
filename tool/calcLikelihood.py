#%%
import matplotlib.pyplot as plt
import numpy as np


def probabilityOfDetectedTumor(x, threshold, sigma):
    # x: true tumor concentration
    # threshold: threshold for detection
    # sigma: width of the sigmoid

    # return: probability of detection with MRI
    return 0.5 + 0.5 * np.sign(x-threshold) * (1. - np.exp(- (x-threshold)**2 / sigma**2))


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

def likelihood(x, xMeasured, stdMeasured):
    #TODO   
    return 1

#dice
def dice(a, b):
    return 2 * np.sum(a * b) / (np.sum(a) + np.sum(b))

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
