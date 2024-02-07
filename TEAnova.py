##TE Anova Setup

import scipy.stats as stats
import numpy as np

def AvgPairWise(data):
    total = 0
    n = 0
    for i in range(len(data)):
        for j in range(i+1,len(data)): 
            total += np.linalg.norm(data[i]-data[j],2)
            n += 1
    return total/n


#select expression protocol
# 0 - Enriched
# 1 - Enhanced
# 2 - Elevated

folder = ["Enriched","Enhanced","Elevated"]
file_name = folder[2]
sample_rate = 10
set_size = 32

##List of Tissues
tissues = ["Adipose", "Adrenal Gland", "Bone Marrow", "Cortex", "Colon", "Duodenum", "Endometrium", "Esophagus", "Fallopian Tube", "Gallbladder", "Heart", "Kidney", "Liver",
         "Lung", "Lymphatic", "Ovary", "Pancreas", "Placenta", "Prostate", "Rectum", "Salivary Gland", "Skeletal Muscle", "Small Intestine", "Smooth Muscle", "Spleen",
         "Stomach", "Testis", "Thyroid", "Tonsil", "Bladder", "Appendix", "Skin"]

##dictionary key = tissue label, value = list of codon bias of genes
expression = {}

with open(file_name + "TEnonStop.csv", "r") as file:
    data = file.readlines()
    for line in data:
        entry = line.strip().split(',')
        tissue = int(entry[0])
        cdb = []
        for i in entry[1:]:
            cdb.append(float(i))
        if tissue in expression.keys():
            expression[tissue].append(np.asarray(cdb))
        else:
            expression[tissue] = []
            expression[tissue].append(np.asarray(cdb))
valid_keys = []

for i in expression.keys():
    if len(expression[i]) > 50:
        valid_keys.append(i)

final = {}

for i in valid_keys:
    data = expression[i]
    sample = []
    for x in range(sample_rate):
        np.random.shuffle(data)
        sample.append(AvgPairWise(data[:set_size]))
    final[tissues[i]] = sample

with open(file_name + "ANOVA.csv", "w") as file:
    print("APD,Group", file = file)
    for tissue in final.keys():
        for sample in final[tissue]:
            print(str(sample) + "," + str(tissue),file = file)
            

stats.f_oneway(*final.values())
    

        
