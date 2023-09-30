#%%
import numpy as np
import matplotlib.pyplot as plt

#%%
path = "/home/jonas/workspace/programs/GliomaSolver/results/2023_09_26-20_07_14_gen_100_loss_bernoulli/results.npy"
dictbern = np.load(path, allow_pickle=True).item()


path = "/home/jonas/workspace/programs/GliomaSolver/results/2023_09_26-20_26_06_gen_100_loss_dice/results.npy"
dictDice = np.load(path, allow_pickle=True).item()



# %%
print(dictbern["opt_params"])
print(dictDice["opt_params"])
# %%
print(dictbern["final_loss"])
print(dictDice["final_loss"])
# %%
print(dictbern["diceT1_75"])
print(dictDice["diceT1_75"])
# %%
print(dictbern["diceFLAIR_30"])
print(dictDice["diceFLAIR_30"])
# %%

def getOpt(path = "/home/jonas/workspace/programs/GliomaSolver/results/2023_09_28-08_46_42_gen_500_loss_bernoulli/results.npy"):

    dict = np.load(path, allow_pickle=True).item()
    #print(dict)
    print('')
    print("params", dict["opt_params"])
    print("final_loss", dict["final_loss"])
    print("diceT1_67", dict["diceT1_67"])
    print("diceFLAIR_25", dict["diceFLAIR_25"])
    #return dict["opt_params"],

    return dict
    
# %%
print('500')
getOpt()
print('1000')
getOpt("/home/jonas/workspace/programs/GliomaSolver/results/2023_09_28-04_08_50_gen_1000_loss_bernoulli/results.npy")
getOpt("/home/jonas/workspace/programs/GliomaSolver/results/2023_09_27-22_16_02_gen_100_loss_bernoulli/results.npy")


# %%
getOpt("/home/jonas/workspace/programs/GliomaSolver/results/2023_09_28-12_00_40_gen_10_loss_bernoulli/results.npy")
# %%
getOpt("/home/jonas/workspace/programs/GliomaSolver/results/2023_09_28-15_45_10_gen_100_loss_dice/results.npy")
# %% this is good! evaluate time... 
dict = getOpt("/home/jonas/workspace/programs/GliomaSolver/results/2023_09_28-14_09_49_gen_100_loss_dice/results.npy")
# %%
getOpt("/home/jonas/workspace/programs/GliomaSolver/results/2023_09_28-16_16_16_gen_50_loss_bernoulli/results.npy")

# %%
getOpt('/home/jonas/workspace/programs/GliomaSolver/resultsSynData/2023_09_28-22_21_27_gen_2_loss_dice/results.npy')
# %%
dict = getOpt("/home/jonas/workspace/programs/GliomaSolver/resultsSynData/patient_0_dtime2023_09_29-03_37_18_gen_100_loss_bernoulli_prior_True_dw0_0.008038116132598106_rho0_0.3065637420701083/results.npy")
# %%
def getDict(path):
    results = np.load(path + '/results.npy', allow_pickle=True).item()

    settings = np.load(path + '/settings.npy', allow_pickle=True).item()

    return results, settings

def plotValues(results, settings):
    import matplotlib.pyplot as plt
    plt.title('values - loss: ' +  settings["lossfunction"] 
              + ' - time: ' + str(round(results["time_min"]/60,1))
    + 'h loss: ' + str(round(results["final_loss"],1)))
    vals = ['x', 'y', 'z', 'D-cm-d', 'rho1-d']
    for i in range(len(vals)):
        plt.plot(results['nsamples'],np.array(results['xmeans']).T[i], label = vals[i])
    plt.yscale('log')
    plt.xlabel('# Samples')
    plt.ylabel('Values')
    plt.legend()
    plt.show()

def plotLoss(results, settings):
    plt.title('loss - loss: ' +  settings["lossfunction"] 
              + ' - time: ' + str(round(results["time_min"]/60,1))
    + 'h loss: ' + str(round(results["final_loss"],1)))

    plt.plot(results['nsamples'],np.array(results['y0s']))
    #plt.yscale('log')
    plt.xlabel('# Samples')
    plt.ylabel('Loss')
    plt.legend()
    plt.show()

def printValues(results, settings, mode = 'abs'):
    vals = ['x', 'y', 'z', 'D-cm-d', 'rho1-d']
    for i in range(len(settings['xMeasured'])):
        if mode == 'rel':
            res = (np.array(results['xmeans']).T[i] - settings['xMeasured'][i] ) /np.array(results['xmeans']).T[i]
            plt.title('relative error')
        else:
            res = (np.array(results['xmeans']).T[i] - settings['xMeasured'][i] )
            plt.title('absolute error')
        plt.plot(results['nsamples'], res, label = vals[i])
        plt.ymin = -1
    plt.ylabel('relative error')
    plt.xlabel('# Samples')
    plt.legend()
    plt.show()

pathes = [#"/home/jonas/workspace/programs/GliomaSolver/resultsSynData/patient_0_dtime2023_09_29-10_09_59_gen_50_loss_dice_prior_True_dw0_0.008038116132598106_rho0_0.3065637420701083", 
          "/home/jonas/workspace/programs/GliomaSolver/resultsSynData/patient_0_dtime2023_09_29-10_03_54_gen_20_loss_dice_prior_True_dw0_0.008038116132598106_rho0_0.3065637420701083", "/home/jonas/workspace/programs/GliomaSolver/resultsSynData/patient_0_dtime2023_09_29-12_01_40_gen_20_loss_bernoulli_prior_True_dw0_0.008038116132598106_rho0_0.3065637420701083/" , "/home/jonas/workspace/programs/GliomaSolver/resultsSynData/patient_0_dtime2023_09_29-13_35_00_gen_10_loss_dice_prior_False_dw0_0.2_rho0_0.1/"]
pathRealPatient = ["/home/jonas/workspace/programs/GliomaSolver/results/2023_09_28-16_16_16_gen_50_loss_bernoulli/", "/home/jonas/workspace/programs/GliomaSolver/results/2023_09_28-14_09_49_gen_100_loss_dice"]

path500 = ["/home/jonas/workspace/programs/GliomaSolver/resultsSynData/patient_0_dtime2023_09_30-06_28_09_gen_500_loss_dice_prior_False_dw0_0.2_rho0_0.1/", "/home/jonas/workspace/programs/GliomaSolver/resultsSynData/patient_0_dtime2023_09_30-03_27_18_gen_500_loss_bernoulli_prior_True_dw0_0.008038116132598106_rho0_0.3065637420701083/", "/home/jonas/workspace/programs/GliomaSolver/resultsSynData/patient_0_dtime2023_09_29-23_44_35_gen_500_loss_dice_prior_True_dw0_0.008038116132598106_rho0_0.3065637420701083/"]

for path in path500:#pathes:pathRealPatient
    results, settings = getDict(path)
    plotLoss(results, settings)
    plotValues(results, settings)
    plt.show()
    print('Final Dice T1', round(results['diceT1_67'],2))
    print('Final Dice FLAIR', round(results['diceFLAIR_25'], 2))
    print('Final Loss', round(results['final_loss'], 2))

    printValues(results, settings)

    #break
    
#plotValues(results, settings)



#%%
results.keys()
results, settings = getDict()
plotLoss(results, settings)
plotValues(results, settings)

# %%
print(dict["loss_function"])
plt.plot(dict['nsamples'], dict['y0s'])# %%
plt.xlabel('nsamples')
plt.ylabel('loss')
dict.keys()
# %%
dict.keys()

# %%
