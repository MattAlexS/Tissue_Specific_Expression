from scipy import stats
from random import sample 

bias = {}


with open("Human101TisGenes.csv", "r") as file:
    data = file.readlines()
    codons = data[0].strip().split(',')[1:]
    for line in data[1:]:
        temp = line.strip().split(',')
        for i in range(len(temp[1:])):
            temp[i+1]=float(temp[i+1])
        bias[temp[0]] = temp[1:]

tissue = {}

with open("TissueExpConcRatio1.5.csv", "r") as file:
    data = file.readlines()
    tis = data[0].strip().split(',')[1:]
    for i in tis:
        tissue[i] = [[],[]]
    for line in data[1:]:
        temp = line.strip().split(',')
        name = temp[0]
        temp = temp[1:]
        for i in range(len(temp)):
            if temp[i] == "1":
                tissue[tis[i]][0].append(name)
            else:
                tissue[tis[i]][1].append(name)
                
for tis in tissue.keys():
    print(tis)
    print(len(tissue[tis][0]))
                
tissues = list(tissue.keys())

n = 300
dataset = {}
for tis in tissue.keys():
    one = sample(tissue[tis][0],n)
    zero = sample(tissue[tis][1],n)
    data1 = []
    data0 = []
    for i in range(len(codons)):
        temp1 = []
        temp0 = []
        for x in range(len(one)):
            temp1.append(bias[one[x]][i])
            temp0.append(bias[zero[x]][i])
        data1.append(temp1)
        data0.append(temp0)
    dataset[tis] = [data1,data0]

p = {}

for tis in dataset.keys():
    row = []
    for codon in range(len(codons)):
        row.append(str(stats.ttest_ind(dataset[tis][0][codon], dataset[tis][1][codon], equal_var = False)[1]))
    p[tis] = row

with open("PC1.5pvaln=" + str(n) + ".csv", "w") as file:
    print("Tissue," + ",".join(codons), file = file)
    for tis in p.keys():
        print(tis + "," + ",".join(p[tis]), file=file)
        
    



        
    
