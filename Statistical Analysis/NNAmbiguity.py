import numpy as np
import pandas as pd

names = ["Adipose", "Adrenal Gland", "Bone Marrow", "Cortex", "Colon", "Duodenum", "Endometrium", "Esophagus", "Fallopian Tube", "Gallbladder", "Heart", "Kidney", "Liver",
         "Lung", "Lymphatic", "Ovary", "Pancreas", "Placenta", "Prostate", "Rectum", "Salivary Gland", "Skeletal Muscle", "Small Intestine", "Smooth Muscle", "Spleen",
         "Stomach", "Testis", "Thyroid", "Tonsil", "Bladder", "Appendix", "Skin"]

included = []

folder = ["Enriched", "Enhanced", "Elevated"]

neighbors = []

neigh_num = 40

protocol = folder[1]

data = pd.read_csv(protocol + "TE.csv", header = None)

with open("NN" + protocol + ".csv", 'w') as file:
    for row in range(len(data)):
        dist_list = [20]*neigh_num
        neigh_list = [0]*neigh_num
        for row2 in range(len(data)):
            if row != row2:
                dist = np.linalg.norm(data.iloc[row,1:] - data.iloc[row2,1:], 2)
                i = 0
                o = 0
                while o==0 and i < neigh_num:
                    if dist < dist_list[i]:
                        dist_list.insert(i,dist)
                        neigh_list.insert(i,str(data.iloc[row2,0]))
                        dist_list.pop()
                        neigh_list.pop()
                        o = 1
                    else:
                        i += 1
        line = (str(data.iloc[row,0]) + "," + ",".join(neigh_list))
        print(line, file = file)


        
            
                
                
            
