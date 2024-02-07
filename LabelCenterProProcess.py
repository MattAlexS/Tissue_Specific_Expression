import pandas as pd
import numpy as np

names = ["Adipose", "Adrenal Gland", "Bone Marrow", "Cortex", "Colon", "Duodenum", "Endometrium", "Esophagus", "Fallopian Tube", "Gallbladder", "Heart", "Kidney", "Liver",
         "Lung", "Lymphatic", "Ovary", "Pancreas", "Placenta", "Prostate", "Rectum", "Salivary Gland", "Skeletal Muscle", "Small Intestine", "Smooth Muscle", "Spleen",
         "Stomach", "Testis", "Thyroid", "Tonsil", "Bladder", "Appendix", "Skin"]

included = []

folder = ["Enriched", "Enhanced", "Elevated"]

poi = {}

protocol = folder[2]

data = pd.read_csv(protocol + "TEnonStop.csv", header = None)

with open(protocol + "CentersNonStop.csv", "w") as file:
    for i in range(32):
        sub = data[data[0] == i].drop(columns = 0)
        if len(sub) > 4:
            included.append(names[i])
            center = sub.mean(axis=0)
            line = [names[i]]
            for i in center:
                line.append(str(i))
            line = ",".join(line)
            print(line, file = file)
        
"""
                        
        poi[names[i]] = [center]
        spread = []
        percent = []
        for row in range(len(sub)):
            spread.append(np.linalg.norm(center - sub.iloc[row,:], 2))
            percent.append(np.percentile(spread,100))
            percent.append(np.percentile(spread,90))
            percent.append(np.percentile(spread,50))
        poi[names[i]].append(percent)

output = []

header = []
header.append("Original Tissue Membership")
for tis in included:
    header.append(tis + "100")
    header.append(tis + "90")
    header.append(tis + "50")

output.append(','.join(header))
    

for row in range(len(data)):
    temp = []
    temp.append(str(names[data.iloc[row,0]]))
    for tissue in included:
        dist = np.linalg.norm(data.iloc[row,1:] - poi[tissue][0], 2)
        if dist <= poi[tissue][1][0]:
            temp.append("1")
        else:
            temp.append("0")

        if dist <= poi[tissue][1][1]:
            temp.append("1")
        else:
            temp.append("0")

        if dist <= poi[tissue][1][2]:
            temp.append("1")
        else:
            temp.append("0")
    output.append(','.join(temp))

with open(protocol + "TisMem.csv", "w") as file:
    for line in output:
        print(line, file = file)
        
"""   
        
        
    
        
