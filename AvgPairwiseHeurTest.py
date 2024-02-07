##Compare Average Pairwise distance to quicker Heuristic

from numpy import random
import statistics
import matplotlib.pyplot as plt
import os
import numpy as np


##Read in Human Codon Bias
#For now ignore duplicates by overwriting
#Revisit appropriate procedure
genes = {}
all_genes = []
with open("Homo_sapiens.GRCh38.codon_bias.csv", "r") as file:
    data = file.readlines()
    for i in data[1:]:
        line = i.strip().split(',')
        name = line[0]
        bias = np.array([float(x) for x in line[1:]])
        all_genes.append(bias)
        genes[name] = bias

def AvgPairWise(data):
    total = 0
    n = 0
    for i in range(len(data)):
        for j in range(i+1,len(data)): 
            total += np.linalg.norm(data[i]-data[j],2)
            n += 1
    return total/n

def HeuristicAPW(data):
    n = len(data)
    z = n
    compares = 0
    total = 0
    sum1 = np.array([0.0]*len(data[0]))

    for point in data:
        for i in range(len(point)):
            sum1[i] += point[i]

    center = list(map(lambda x : x/n, sum1))
    data = sorted(data, key = lambda x: np.linalg.norm(x - center,2), reverse = True)
    for point in data[:-1]:
        sum1 -= point
        z -= 1
        compares += z
        center = np.true_divide(sum1, z)
        total += (np.linalg.norm(point - center, 2) * z)
    return total/compares
    

full_size = len(all_genes)

num_trials = 50
set_size = 50
avgs = []
heurs = []

for x in range(num_trials):
    sample = []
    for i in range(set_size):
        sample.append(all_genes[random.randint(0,full_size)])
    heurs.append(HeuristicAPW(sample))
    avgs.append(AvgPairWise(sample))

x = range(num_trials)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(x, avgs, s=10, c='b', marker="s", label='Average Pairwise')
ax1.scatter(x, heurs, s=10, c='r', marker="o", label='Heuristic')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()
    


