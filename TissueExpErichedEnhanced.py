##Inestigate implication of codon usage on tissue expression
##Compare codon bias to tissue specific tRNA concentrations
##Use Tissue expression data to determine which genes are expressed in which tissues

##Processing input files
##Hash for converting Ensembl Identifies to UniProtKB
import statistics
import matplotlib.pyplot as plt
import os
import numpy as np

key = {}

with open("Ensembl-UniProt.csv", "r") as file:
    data = file.readlines()
    for i in data:
        line = i.strip().split(",")
        key[line[0]] = [line[1],line[2]]

##Reading in Tissue Expression Input       
entries = []
with open("TissueExp.tsv", 'r') as file:
    data = file.readlines()
    for i in range(len(data)):
        if data[i][0] != "#":
            temp = data[i].strip().split("\t")
            for x in range(len(temp)):
                if temp[x] == '':
                    temp[x] = 0
                elif i > 4 and x > 1:
                    temp[x] = float(temp[x])
            extra = 34 - len(temp)
            for addition in range(extra):
                temp.append(0)
            entries.append(temp)

"""
entries = []

with open("TissueExp.tsv", 'r') as file:
    data = file.readlines()
    for i in data:
        if i[0] != "#":
            temp = i.strip().split("\t")
            for i in range(len(temp)):
                if temp[i] == '':
                    temp[i] = 0
            extra = 34 - len(temp)
            for i in range(extra):
                temp.append(0)
            entries.append(temp)
"""

##Sorting Genes according to Tissue they are most highly expressed in

names = ["Adipose", "Adrenal Gland", "Bone Marrow", "Cortex", "Colon", "Duodenum", "Endometrium", "Esophagus", "Fallopian Tube", "Gallbladder", "Heart", "Kidney", "Liver",
         "Lung", "Lymphatic", "Ovary", "Pancreas", "Placenta", "Prostate", "Rectum", "Salivary Gland", "Skeletal Muscle", "Small Intestine", "Smooth Muscle", "Spleen",
         "Stomach", "Testis", "Thyroid", "Tonsil", "Bladder", "Appendix", "Skin"]
human_enr = []
human_enh = []
human_ele = []

for i in names:
    human_enr.append([])
    human_enh.append([])
    human_ele.append([])

for i in range(1,len(entries)):
    if entries[i][0] in key.keys() and key[entries[i][0]][0] != '':
        idtag = key[entries[i][0]][0]
        gene_mean = statistics.mean(entries[i][2:])
        enhanced = gene_mean * 5
        temp = entries[i][2:]
        top = max(temp)
        temp.remove(top)
        enriched = max(temp) * 5

        for x in range(len(human_enh)):
            if float(entries[i][x+2]) > enhanced:
                human_enh[x].append(idtag)
            if float(entries[i][x+2]) > enriched:
                human_enr[x].append(idtag)
            if float(entries[i][x+2]) == top:
                human_ele[x].append(idtag)
"""
for i in range(1,len(entries)):
    if entries[i][0] in key.keys() and key[entries[i][0]][0] != '':
        idtag = key[entries[i][0]][0]
        gene_mean = statistics.mean(entries[i][2:])
        gene_std = statistics.stdev(entries[i][2:])
        five_sigma = 5 * gene_std
        enhanced = gene_mean + five_sigma
        for x in range(len(human_enh)):
            if float(entries[i][x+2]) > enhanced:
                human_enh[x].append(idtag)
            expression = float(entries[i][x+2]) - five_sigma
            enriched = True
            for t in range(len(human_enr)):
                if expression < float(entries[i][t+2]):
                    enriched = False
            if enriched == True:
                human_enr[x].append(idtag)
        human[tissue].append(idtag)
"""

##Distance Function
def EuclidianDistance(point,cluster):
    if len(point) != len(cluster):
        return False
    else:
        diff = []
        for i in range(len(point)):
            point[i] = float(point[i])
            cluster[i] = float(cluster[i])                            
            diff.append((point[i] - cluster[i])**2)
        dist = (sum(diff))**0.5
        return dist


##Read in Codon data for human genes
##Compute averages for multiple copies of same gene
genes = {}
temp = []
with open("Homo_sapiens.GRCh38.codon_bias.csv", "r") as file:
    data = file.readlines()
    for i in data[1:]:
        line = i.strip().split(',')
        name = line[0]
        temp = []
        for x in line[1:]:
            temp.append(float(x))
        genes[name] = np.asarray(temp)
"""
genes = {}
temp = []
with open("Homo_sapiens.GRCh38.codon_bias.csv", "r") as file:
    data = file.readlines()
    for i in data[1:]:
        line = i.strip().split(',')
        name = line[0]
        bias = line[1:]
        if name in genes.keys():
            new = []
            multiple = genes[name][0]
            for i in range(len(genes[name][1])):
                new.append(((float(name[1][i])*multiple) + float(bias[i]))/(multiple + 1))
            genes[name][0] += 1
            genes[name][1] = new
        else:
            genes[name] = [1,bias]
"""
"""       
for i in temp:
    name = i[0]
    bias = i[1:]
    collection = []
    collection.append(bias)
    temp.remove(i)
    for x in temp:
        if x[0] == name:
            collection.append(x[1:])
            temp.remove(x)
    if len(collection) > 1:
        average = []
        for i in range(len(collection[0])):
            total = 0
            for x in range(len(collection)):
                total += float(collection[x][i])
            average.append(total/len(collection))
        genes[name] = average
    else:
        genes[name] = bias

##Directory Setup
script_dir = os.path.dirname(__file__)
results_dir = os.path.join(script_dir, 'Tissue Elevated Results/')
"""

##Calculate pairwise distances
summary = []
for i in range(len(human_enr)):
    print(names[i])
    total = str(len(human_enr[i]))
    sumline = [names[i],str(2.114690180703266),str(0.5753425783877906),str(len(human_enr[i]))]
    pair_dist_total = 0
    compares = 0
    for j in range(len(human_enr[i]) - 1, -1, -1):
        for k in range(j):
            if human_enr[i][j] in genes.keys() and human_enr[i][k] in genes.keys():
                pair_dist_total += np.linalg.norm(genes[human_enr[i][j]] - genes[human_enr[i][k]],2)
                compares += 1
    if compares == 0:
        sumline.append(str(0))
    else:
        sumline.append(str(pair_dist_total/compares))
    summary.append(sumline)

with open("HumanEnrichedSignificanceTest.csv", "w") as file:
    for i in summary:
        print(",".join(i), file = file)

    
                                                
            
            
    
         
                









    
