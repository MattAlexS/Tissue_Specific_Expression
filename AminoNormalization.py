import statistics
import matplotlib.pyplot as plt
import os
import numpy as np
import random


aminos = {
    'Lys':[0,2],
    'Asn':[1,3],
    'Thr':[4,5,6,7],
    'Arg':[8,10,23,24,25,26],
    'Ser':[9,11,51,52,53,54],
    'Ile':[12,13,14],
    'Gln':[15,17],
    'His':[16,18],
    'Pro':[19,20,21,22],
    'Leu':[27,28,29,30,58,60],
    'Glu':[31,33],
    'Asp':[32,34],
    'Ala':[35,36,37,38],
    'Gly':[39,40,41,42],
    'Val':[43,44,45,46],
    'Stop':[47,49,55],
    'Tyr':[48,50],
    'Cys':[56,57],
    'Phe':[59,61]
    }

folder = ["Enriched", "Enhanced", "Elevated"]
protocol = folder[0]
metric_num = 30
sample_rate = 50

labels = []
genes = []
with open(protocol + "TE.csv", "r") as file:
    data = file.readlines()
    for i in data:
        line = i.strip().split(',')
        temp = []
        labels.append(line[0])
        for x in line[1:]:
            temp.append(float(x))
        genes.append(temp)

two_labels = []
two_apd =[]

three_labels = []
three_apd = []

four_labels = []
four_apd = []

six_labels = []
six_apd = []


for acid in aminos.keys():
    meanlist = []
    for x in range(sample_rate):
        sample = random.sample(genes,k=metric_num)
        subset = []
        for unit in sample:
            temp = []
            for index in aminos[acid]:
                temp.append(unit[index])
            subset.append(np.asarray(temp))
        pair_dist_total = 0
        compares = 0
        for j in range(metric_num -1, -1, -1):
            for k in range(j):
                pair_dist_total += np.linalg.norm(subset[j] - subset[k],2)
                compares += 1
        meanlist.append(pair_dist_total/compares)
    if len(aminos[acid]) == 2:
        if len(two_apd) == 0:
            two_apd.append(meanlist)
            two_labels.append(acid)
        else:
            rank = 0
            while rank < len(two_apd):
                if statistics.mean(meanlist) > statistics.mean(two_apd[rank]):
                    rank += 1
                else:
                    break
            two_apd.insert(rank, meanlist)
            two_labels.insert(rank, acid)
    elif len(aminos[acid]) == 3:
        if len(three_apd) == 0:
            three_apd.append(meanlist)
            three_labels.append(acid)
        else:
            rank = 0
            while rank < len(three_apd):
                if statistics.mean(meanlist) > statistics.mean(three_apd[rank]):
                    rank += 1
                else:
                    break
            three_apd.insert(rank, meanlist)
            three_labels.insert(rank, acid)
    elif len(aminos[acid]) == 4:
        if len(four_apd) == 0:
            four_apd.append(meanlist)
            four_labels.append(acid)
        else:
            rank = 0
            while rank < len(four_apd):
                if statistics.mean(meanlist) > statistics.mean(four_apd[rank]):
                    rank += 1
                else:
                    break
            four_apd.insert(rank, meanlist)
            four_labels.insert(rank, acid)
    elif len(aminos[acid]) == 6:
        if len(six_apd) == 0:
            six_apd.append(meanlist)
            six_labels.append(acid)
        else:
            rank = 0
            while rank < len(six_apd):
                if statistics.mean(meanlist) > statistics.mean(six_apd[rank]):
                    rank += 1
                else:
                    break
            six_apd.insert(rank, meanlist)
            six_labels.insert(rank, acid)

meanAPD = {}

maxAPD = 0

for i in range(len(two_apd)):
    meanAPD[two_labels[i]] = statistics.mean(two_apd[i])
    if statistics.mean(two_apd[i]) > maxAPD:
        maxAPD = statistics.mean(two_apd[i])
"""
for i in range(len(three_apd)):
    meanAPD[three_labels[i]] = statistics.mean(three_apd[i])
    if statistics.mean(three_apd[i]) > maxAPD:
        maxAPD = statistics.mean(three_apd[i])
"""
meanAPD[three_labels[0]] = statistics.mean(three_apd[0])
if statistics.mean(three_apd[0]) > maxAPD:
    maxAPD = statistics.mean(three_apd[0])

for i in range(len(four_apd)):
    meanAPD[four_labels[i]] = statistics.mean(four_apd[i])
    if statistics.mean(four_apd[i]) > maxAPD:
        maxAPD = statistics.mean(four_apd[i])

for i in range(len(six_apd)):
    meanAPD[six_labels[i]] = statistics.mean(six_apd[i])
    if statistics.mean(six_apd[i]) > maxAPD:
        maxAPD = statistics.mean(six_apd[i])

weights = {}

for amino in meanAPD.keys():
    weights[amino] = maxAPD/meanAPD[amino]

weights["Stop"] = 0.0

print(weights)

with open(protocol + "NormalizationTE.csv","w") as file:
    for i in range(len(genes)):
        for amino in aminos.keys():
            for index in aminos[amino]:
                genes[i][index] = str(genes[i][index] * weights[amino])
        print(str(labels[i]) + "," + ",".join(genes[i]),file = file)
        
    

"""        
plt.boxplot(three_apd, labels = three_labels)
plt.title("Average Pairwise Distance Within Amino Acids with 3 Codons (" + protocol + " Expression)\nSample Rate = " + str(sample_rate))
plt.ylabel("Average Pairwise Distance (Sample Rate = " + str(metric_num) + ")")
plt.xlabel("Amino Acid")
#plt.axhline(y=background, color='r', linestyle='-')
plt.xticks(rotation=90)
plt.show()   
"""




    
