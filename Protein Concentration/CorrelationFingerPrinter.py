import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import random

matplotlib.style.use('ggplot')

names = ['Adrenal Gland', 'Cerebral Cortex', 'Colon', 'Esophagus', 'Gallbladder', 'Heart', 'Kidney',
         'Liver', 'Lung', 'Ovary', 'Pancreas', 'Placenta', 'Prostate', 'Rectum', ' Testis', 'Urinary Bladder',"None"]


acids = {
    'Lys':[1,3],
    'Asn':[2,4],
    'Thr':[5,6,7,8],
    'Arg':[9,11,24,25,26,27],
    'Ser':[10,12,52,53,54,55],
    'Ile':[13,14,15],
    'Gln':[16,18],
    'His':[17,19],
    'Pro':[20,21,22,23],
    'Leu':[28,29,30,31,59,61],
    'Glu':[32,34],
    'Asp':[33,35],
    'Ala':[36,37,38,39],
    'Gly':[40,41,42,43],
    'Val':[44,45,46,47],
    'Stop':[48,50,56],
    'Tyr':[49,51],
    'Cys':[57,58],
    'Phe':[60,62]
    }

percentages = []
"""
folder = ["Enriched", "Enhanced", "Elevated"]

protocol = folder[1]
"""
data = pd.read_csv("TisTRatio1.5Classes.csv", header = None)

#train:test ratio
ratio = 6
featnum = 4

model = {}

test = data.iloc[range(0,len(data),ratio),:]

train = data.drop(index = range(0,len(data),ratio))

for i in range(len(names)):
    sub = train[train[0] == i].drop(columns = 0)
    #if len(sub) > 5:
    model[names[i]] = zeros = [[0 for a in range(len(sub.columns))] for b in range(len(sub.columns))]
    for j in range(len(sub.columns)-1):
        for k in range(j+1,len(sub.columns)):
            X = pd.DataFrame(sub.iloc[:,j])
            y = pd.DataFrame(sub.iloc[:,k])
            reg = LinearRegression().fit(X,y)
            model[names[i]][j][k] = [reg, reg.score(X,y)]

    #write test set routine
for z in range(10):
    cut = random.sample(range(1,62),(62 - featnum))
    invcut = []
    for c in range(1,62):
        if c not in cut:
            invcut.append(c)
            
    feat = test.drop(columns = cut)
    #feat.columns = range(feat.shape[1])

    output = []
    
    for i in range(len(feat)):
        label = names[feat.iloc[i,0]]
        example = feat.iloc[i,1:]
        maximum = 0
        predict = False
        for tissue in model.keys():
            vote = 0
            for j in range(1,len(invcut)-1):
                for k in range(j+1,len(invcut)):
                    guess = model[tissue][invcut[j]][invcut[k]][0].predict(example[invcut[j]].reshape(1,-1))
                    actual = example[invcut[k]]
                    corr =  model[tissue][invcut[j]][invcut[k]][1]
                    if guess > 0.5:
                        maxDiff = guess
                    else:
                        maxDiff = 1 - guess
                    vote += ((maxDiff - (abs(guess - actual)))/maxDiff)*corr
            if vote > maximum:
                maximum = vote
                predict = tissue
        output.append([label,predict])
    total = 0
    for p in output:
        if p[0] == p[1]:
            total += 1
    percentages.append((total/len(output))*100)
    
    
    
            
                    
                
