##Distribution of pairwise distances within a point set of varying size.
##Point Set size 2 --> 2048
##Measure Mean and Variance of pairwise distances
##Data Set: Human Codon Bias

##Libraries

from numpy import random
import statistics
import matplotlib.pyplot as plt
import os
import numpy as np



##Read in Human Codon Bias
#For now ignore duplicates by overwriting
#Revisit appropriate procedure
genes = {}
with open("Homo_sapiens.GRCh38.codon_bias.csv", "r") as file:
    data = file.readlines()
    for i in data[1:]:
        line = i.strip().split(',')
        name = line[0]
        bias = line[1:]
        genes[name] = bias

data = []
for i in genes.keys():
    temp = []
    for col in genes[i]:
        temp.append(float(col))
    data.append(np.asarray(temp))
    
"""
script_dir = os.path.dirname(__file__)
results_dir = os.path.join(script_dir, 'Gene Set Size Study/')

for power in range(2, 12):
    set_size = 2**power
    means = []
    for i in range(30):
        sample = []
        dists=[]
        sample.append(data[random.randint(0,len(data))])
        for x in range(set_size - 1):
            pick = data[random.randint(0,len(data))]
            for entry in sample:
                dists.append(EuclidianDistance(pick,entry))
            sample.append(pick)
        means.append(statistics.mean(dists))
    gmean = statistics.mean(means)
    stdev = statistics.stdev(means)
    plt.hist(means)
    plt.title("Mean Pairwise Distance at Set Size: " + str(set_size) + "\nMean: " + str(gmean) + "\nStandard Deviation: " + str(stdev))
    plt.xlabel("Mean Pairwise Distances")
    plt.ylabel("Frequency")
    plt.savefig(fname = results_dir + str(set_size), dpi = 900, bbox_inches = "tight", pad_inches = 0.5)
    plt.close()
"""
print("Assembled")


def AvgPairWise(data):
    total = 0
    n = 0
    for i in range(len(data)):
        for j in range(i+1,len(data)): 
            total += np.linalg.norm(data[i]-data[j],2)
            n += 1
    return total/n



"""
mu = total/n
ex2 = sqtotal/n

var = sqtotal - mu**2

mu = 2.114690180703266

sd = 0.5753425783877906
"""


    
            
            
