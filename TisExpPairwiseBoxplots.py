##Inestigate implication of codon usage on tissue expression
##Compare codon bias to tissue specific tRNA concentrations
##Use Tissue expression data to determine which genes are expressed in which tissues

##Processing input files
##Hash for converting Ensembl Identifies to UniProtKB
import statistics
import matplotlib.pyplot as plt
import os
import numpy as np
import random
import datetime

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
        for x in range(len(line)):
            if x == 0 or x == 48 or x == 50 or x == 56:
                continue
            else:
                temp.append(float(line[x]))
        genes[name] = np.asarray(temp)




full_list = []
for i in genes.keys():
    full_list.append(i)
background_total = 0

study = human_ele
metric_num = 30
sample_rate = 1000
back = []
labels = []
summary = []

for i in range(sample_rate):
    pair_dist_total = 0
    compares = 0
    np.random.shuffle(full_list)
    sample = full_list[:metric_num]
    for j in range(metric_num -1, -1, -1):
        for k in range(j):
            pair_dist_total += np.linalg.norm(genes[sample[j]] - genes[sample[k]],2)
            compares += 1
    background_total += (pair_dist_total/compares)
    back.append(pair_dist_total/compares)

background = background_total/sample_rate

summary.append(back)
labels.append("Background")


##Calculate pairwise distances

for i in range(len(study)):
    print(names[i])
    for j in range(len(study[i])-1, -1, -1):
        if study[i][j] not in genes.keys():
            del study[i][j]
    if len(study[i]) > metric_num:
        meanlist = []
        for x in range(sample_rate):
            pair_dist_total = 0
            compares = 0
            np.random.shuffle(study[i])
            sample = random.sample(study[i],k=metric_num)
            for j in range(metric_num -1, -1, -1):
                for k in range(j):
                    pair_dist_total += np.linalg.norm(genes[sample[j]] - genes[sample[k]],2)
                    compares += 1
            meanlist.append(pair_dist_total/compares)
        if len(summary) == 0:
            summary.append(meanlist)
            labels.append(names[i])
        else:
            rank = 0
            while rank < len(summary):
                if statistics.mean(meanlist) > statistics.mean(summary[rank]):
                    rank += 1
                else:
                    break
            summary.insert(rank, meanlist)
            labels.insert(rank, names[i])



        


plt.boxplot(summary, labels = labels)
plt.title("Average Pairwise Distance Within Tissue (Elevated Expression)\nSample Rate = " + str(sample_rate))
plt.ylabel("Average Pairwise Distance (Sample Rate = " + str(metric_num) + ")")
plt.xlabel("Tissue Type")
plt.axhline(y=background, color='r', linestyle='-')
plt.xticks(rotation=90)
plt.show()

"""
with open("ElevatedTE.csv","w") as file:
    for i in range(len(human_ele)):
        for gene in human_ele[i]:
            if gene in genes.keys():
                new = genes[gene].tolist()
                new.insert(0,str(i))
                for j in range(len(new)):
                    new[j]=str(new[j])
                print(",".join(new), file = file)
"""
                                                
            
            
    
         
                









    
