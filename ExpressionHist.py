#Tissue Expression Within Tissue

import numpy as np
import matplotlib.pyplot as plt

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
            
adipose = []
adrenal = []
marrow = []
cortex = []
colon = []
duodenum = []
endometrium = []
esophagus = []
fallopian = []
gallbladder = []
heart = []
kidney = []
liver = []
lung = []
lymph = []
ovary = []
pancreas = []
placenta = []
prostate = []
rectum = []
salivary = []
skeletalmuscle = []
smallintestine = []
smoothmuscle = [] 
spleen = []
stomach = []
testis = []
thyroid = []
tonsil = []
bladder = []
appendix = []
skin = []

misses = 0

for i in entries[1:]:
    if len(i) == 34:
        adipose.append(float(i[2]))
        adrenal.append(float(i[3]))
        marrow.append(float(i[4]))
        cortex.append(float(i[5]))
        colon.append(float(i[6]))
        duodenum.append(float(i[7]))
        endometrium.append(float(i[8]))
        esophagus.append(float(i[9]))
        fallopian.append(float(i[10]))
        gallbladder.append(float(i[11]))
        heart.append(float(i[12]))
        kidney.append(float(i[13]))
        liver.append(float(i[14]))
        lung.append(float(i[15]))
        lymph.append(float(i[16]))
        ovary.append(float(i[17]))
        pancreas.append(float(i[18]))
        placenta.append(float(i[19]))
        prostate.append(float(i[20]))
        rectum.append(float(i[21]))
        salivary.append(float(i[22]))
        skeletalmuscle.append(float(i[23]))
        smallintestine.append(float(i[24]))
        smoothmuscle.append(float(i[25]))
        spleen.append(float(i[26]))
        stomach.append(float(i[27]))
        testis.append(float(i[28]))
        thyroid.append(float(i[29]))
        tonsil.append(float(i[30]))
        bladder.append(float(i[31]))
        appendix.append(float(i[32]))
        skin.append(float(i[33]))
    else:
        misses += 1

human = [adipose, adrenal, marrow, cortex, colon, duodenum, endometrium, esophagus, fallopian, gallbladder, heart, kidney, liver,
         lung, lymph, ovary, pancreas, placenta, prostate, rectum, salivary, skeletalmuscle, smallintestine, smoothmuscle, spleen,
         stomach, testis, thyroid, tonsil, bladder, appendix, skin]
names = ["adipose", "adrenal", "marrow", "cortex", "colon", "duodenum", "endometrium", "esophagus", "fallopian", "gallbladder", "heart", "kidney", "liver",
         "lung", "lymph", "ovary", "pancreas", "placenta", "prostate", "rectum", "salivary", "skeletalmuscle", "smallintestine", "smoothmuscle", "spleen",
         "stomach", "testis", "thyroid", "tonsil", "bladder", "appendix", "skin"]

plt.hist(liver, bins = 500)
plt.ylabel("Frequency")
plt.xlabel("TPM")
plt.show()
