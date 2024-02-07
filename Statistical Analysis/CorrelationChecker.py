import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.style.use('ggplot')

names = ["Adipose", "Adrenal Gland", "Bone Marrow", "Cortex", "Colon", "Duodenum", "Endometrium", "Esophagus", "Fallopian Tube", "Gallbladder", "Heart", "Kidney", "Liver",
         "Lung", "Lymphatic", "Ovary", "Pancreas", "Placenta", "Prostate", "Rectum", "Salivary Gland", "Skeletal Muscle", "Small Intestine", "Smooth Muscle", "Spleen",
         "Stomach", "Testis", "Thyroid", "Tonsil", "Bladder", "Appendix", "Skin"]

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


folder = ["Enriched", "Enhanced", "Elevated"]

protocol = folder[1]

data = pd.read_csv(protocol + "TE.csv", header = None)

for i in range(32):
    sub = data[data[0] == i].drop(columns = 0)
    if len(sub) > 5:
        plt.matshow(sub.corr())
        plt.xticks(range(len(sub.columns)), sub.columns)
        plt.yticks(range(len(sub.columns)), sub.columns)
        plt.colorbar()
        plt.show()

"""
        for amino1 in acids.keys():
            for amino2 in acids.keys():
                if amino1 =! amino2:
                    for j in acids

"""
