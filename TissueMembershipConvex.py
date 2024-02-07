import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay


def in_hull(p, hull):
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)
    return hull.find_simplex(p)>=0


names = ["Adipose", "Adrenal Gland", "Bone Marrow", "Cortex", "Colon", "Duodenum", "Endometrium", "Esophagus", "Fallopian Tube", "Gallbladder", "Heart", "Kidney", "Liver",
         "Lung", "Lymphatic", "Ovary", "Pancreas", "Placenta", "Prostate", "Rectum", "Salivary Gland", "Skeletal Muscle", "Small Intestine", "Smooth Muscle", "Spleen",
         "Stomach", "Testis", "Thyroid", "Tonsil", "Bladder", "Appendix", "Skin"]

acids = {
    #'Lys':[1,3],
    #'Asn':[2,4],
    'Thr':[5,6,7],#,8],
    'Arg':[9,11,24,25,26],#,27],
    'Ser':[10,12,52,53,54],#,55],
    'Ile':[13,14],#,15],
    #'Gln':[16,18],
    #'His':[17,19],
    'Pro':[20,21,22],#,23],
    'Leu':[28,29,30,31,59],#,61],
    #'Glu':[32,34],
    #'Asp':[33,35],
    'Ala':[36,37,38],#,39],
    'Gly':[40,41,42,],#43],
    'Val':[44,45,46,],#47],
    'Stop':[48,50]#,56]#,
    #'Tyr':[49,51],
    #'Cys':[57,58],
    #'Phe':[60,62]
    }

included = []

folder = ["Enriched", "Enhanced", "Elevated"]

poi = {}

protocol = folder[2]

data = pd.read_csv(protocol + "TE.csv", header = None)

for i in range(32):
    sub = data[data[0] == i].drop(columns = 0)
    if len(sub) > 5:
        simplex = []
        included.append(names[i])
        for amino in acids:
            hplane = sub[acids[amino]]
            check = data[acids[amino]]
            simplex.append(in_hull(check, hplane))
        poi[names[i]] = simplex

output = []

header = []
header.append("Original Tissue Membership")
for tis in included:
    header.append(tis)

output.append(','.join(header))
    

for row in range(len(data)):
    temp = []
    temp.append(str(names[data.iloc[row,0]]))
    for tissue in included:
        aa = []
        for mem in poi[tissue]:
            aa.append(mem[row])
        if False in aa:
            temp.append("0")
        else:
            temp.append("1")
    output.append(','.join(temp))

with open(protocol + "TisMemCvexHull.csv", "w") as file:
    for line in output:
        print(line, file = file)

            
        
        
    
        
