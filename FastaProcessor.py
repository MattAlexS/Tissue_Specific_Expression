##Process raw Fasta File to extract required sequences

from numpy import random
import statistics

def codonBiasVectorize(seq):
    alphabet = ["A", "C", "G", "T"]
    lencheck = len(seq) % 3
    if lencheck != 0:
        extend = (3 - lencheck)
        seq = seq + "X"*extend
    unknown = 0
    aminos = {
        'K':[0,2],
        'N':[0,2],
        'T':[0,4],
        'R':[0,6],
        'S':[0,6],
        'I':[0,3],
        'M':[0,1],
        'Q':[0,2],
        'H':[0,2],
        'P':[0,4],
        'L':[0,6],
        'E':[0,2],
        'D':[0,2],
        'A':[0,2],
        'G':[0,4],
        'V':[0,4],
        'Stop':[0,3],
        'Y':[0,2],
        'C':[0,2],
        'W':[0,1],
        'F':[0,2]
    }
            
    codons = {
        'AAA':[0,'K'],    #K
        'AAC':[0,'N'],    #N
        'AAG':[0,'K'],    #K
        'AAT':[0,'N'],    #N
        'ACA':[0,'T'],    #T
        'ACC':[0,'T'],    #T
        'ACG':[0,'T'],    #T
        'ACT':[0,'T'],    #T
        'AGA':[0,'R'],    #R
        'AGC':[0,'S'],    #S
        'AGG':[0,'R'],    #R
        'AGT':[0,'S'],    #S
        'ATA':[0,'I'],    #I
        'ATC':[0,'I'],    #I
        'ATG':[0,'M'],    #M
        'ATT':[0,'I'],    #I
        'CAA':[0,'Q'],    #Q
        'CAC':[0,'H'],    #H
        'CAG':[0,'Q'],    #Q
        'CAT':[0,'H'],    #H
        'CCA':[0,'P'],    #P
        'CCC':[0,'P'],    #P
        'CCG':[0,'P'],    #P
        'CCT':[0,'P'],    #P
        'CGA':[0,'R'],    #R
        'CGC':[0,'R'],    #R
        'CGG':[0,'R'],    #R
        'CGT':[0,'R'],    #R
        'CTA':[0,'L'],    #L
        'CTC':[0,'L'],    #L
        'CTG':[0,'L'],    #L
        'CTT':[0,'L'],    #L
        'GAA':[0,'E'],    #E
        'GAC':[0,'D'],    #D
        'GAG':[0,'E'],    #E
        'GAT':[0,'D'],    #D
        'GCA':[0,'A'],    #A
        'GCC':[0,'A'],    #A
        'GCG':[0,'A'],    #A
        'GCT':[0,'A'],    #A
        'GGA':[0,'G'],    #G
        'GGC':[0,'G'],    #G
        'GGG':[0,'G'],    #G
        'GGT':[0,'G'],    #G
        'GTA':[0,'V'],    #V
        'GTC':[0,'V'],    #V
        'GTG':[0,'V'],    #V
        'GTT':[0,'V'],    #V
        'TAA':[0,'Stop'],    #Stop
        'TAC':[0,'Y'],    #Y
        'TAG':[0,'Stop'],    #Stop
        'TAT':[0,'Y'],    #Y
        'TCA':[0,'S'],    #S
        'TCC':[0,'S'],    #S
        'TCG':[0,'S'],    #S
        'TCT':[0,'S'],    #S
        'TGA':[0,'Stop'],    #Stop
        'TGC':[0,'C'],    #C
        'TGG':[0,'W'],    #W
        'TGT':[0,'C'],    #C
        'TTA':[0,'L'],    #L
        'TTC':[0,'F'],    #F
        'TTG':[0,'L'],    #L
        'TTT':[0,'F']    #F
        }
    cursor = 0
    while cursor <= len(seq)-2:
        if seq[cursor] in alphabet and seq[cursor + 1] in alphabet and seq[cursor + 2] in alphabet:
            codons[seq[cursor:cursor + 3]][0] += 1
            aminos[codons[seq[cursor:cursor + 3]][1]][0] += 1
            cursor += 3
        else:
            unknown += 1
            cursor += 3
    codonbias = []
    for i in codons.keys():
        if aminos[codons[i][1]][1] != 1:
            if aminos[codons[i][1]][0] == 0:
                codonbias.append(str(1/aminos[codons[i][1]][1]))
            else:
                codonbias.append(str(codons[i][0]/aminos[codons[i][1]][0]))
    return codonbias

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
######################################################################
fasta = []

with open("Homo_sapiens.GRCh38.cds.all.fa", "r") as file:
    data = file.readlines()
    for i in data:
        line = i.strip()
        if line[0] == ">":
            fasta.append(line.split(' ')[0] + '|' + line.split(' ')[3].split(':')[1] + ">")
        else:
            fasta.append(line)

temp = ''.join(fasta)

fasta = temp.split('>')
fasta.pop(0)
del data
del temp


genes = {}

for i in range(len(fasta)):
    if fasta[i][0] == "E":
        geneID = fasta[i].split("|")[1].split(".")[0]
        transcriptID = fasta[i].split("|")[0].split(".")[0]
        sequence = fasta[i+1]
        if geneID in genes.keys():
            if sequence not in genes[geneID][1]:
                genes[geneID][0].append(transcriptID)
                genes[geneID][1].append(sequence)
        else:
            genes[geneID] = [[transcriptID],[sequence]]
            
del fasta

key = {}

with open("Ensembl-UniProt.csv", "r") as file:
    data = file.readlines()
    for i in data:
        line = i.strip().split(",")
        key[line[0]] = [line[1],line[2]]

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

            
##Sorting Genes according to Tissue Enriched and Tissue Enhanced
names = ["adipose", "adrenal", "marrow", "cortex", "colon", "duodenum", "endometrium", "esophagus", "fallopian", "gallbladder", "heart", "kidney", "liver",
         "lung", "lymph", "ovary", "pancreas", "placenta", "prostate", "rectum", "salivary", "skeletalmuscle", "smallintestine", "smoothmuscle", "spleen",
         "stomach", "testis", "thyroid", "tonsil", "bladder", "appendix", "skin"]
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
                


genesCodonBias = {}

for tissue in human:
    for gene in tissue:
        vectors = []
        for i in genes[gene][1]:
            vectors.append(codonBiasVectorize(i))
        genesCodonBias[gene] = vectors

##Averaging alternate transcripts of the same gene
"""
genesAvg = {}

for i in genesCodonBias.keys():
    num = len(genesCodonBias[i])
    avg = []
    for x in range(len(genesCodonBias[i][0])):
        total = 0
        for y in genesCodonBias[i]:
            total += float(y[x])
        avg.append(total/num)
    genesAvg[i] = avg
"""

##Using the closest Transcript to the moving average of the Tissue
from operator import itemgetter
genesClosest = {}


for tissue in human:
    temp = []
    for gene in tissue:
        temp.append((gene, len(genesCodonBias[gene])))
    ordered = sorted(temp, key=itemgetter(1))
    rollingAvg = []
    for index in range(len(ordered)):
        if index == 0:
            if ordered[index][1] == 1:
                rollingAvg = genesCodonBias[ordered[index][0]][0]
                genesClosest[ordered[index][0]] = rollingAvg
            else:
                for i in range(len(genesCodonBias[ordered[index][0]][0])):
                    codonSum = 0
                    for x in genesCodonBias[ordered[index][0]]:
                        codonSum += float(x[i])
                    rollingAvg.append(codonSum/len(genesCodonBias[ordered[index][0]]))
                genesClosest[ordered[index][0]] = rollingAvg
        else:
            closest = 0
            dist = 400
            for transcript in range(len(genesCodonBias[ordered[index][0]])):
                if EuclidianDistance(genesCodonBias[ordered[index][0]][transcript],rollingAvg) < dist:
                    newdist = EuclidianDistance(genesCodonBias[ordered[index][0]][transcript],rollingAvg)
                    dist = newdist
                    closest = transcript
            genesClosest[ordered[index][0]] = genesCodonBias[ordered[index][0]][closest]
            newAvg = []
            for i in range(len(rollingAvg)):
                codonSum = 0
                codonSum += float(genesClosest[ordered[index][0]][i])
                codonSum += rollingAvg[i] * (index)
                newAvg.append(codonSum/(index+1))
            rollingAvg = newAvg

            

##Generating the Center point of each Tissue
centers = []

for tissue in human:
    temp = []
    for gene in tissue:
        temp.append(genesClosest[gene])
    num = len(temp)
    center = []
    for i in range(len(temp[0])):
        total = 0
        for x in temp:
            total += x[i]
        center.append(total/num)
    centers.append(center)

##Calculating average gene distance to tissue center
averageDist = []

for i in range(len(human)):
    temp = []
    for gene in human[i]:
        temp.append(EuclidianDistance(genesClosest[gene],centers[i]))
    averageDist.append(sum(temp)/len(temp))


##Calculating likelihood of finding and Average Distance lower or as low
allGenesList = []
for i in genesClosest.keys():
    allGenesList.append(i)

pval = []
n = 1000
for i in range(len(human)):
    closer = 0
    for x in range(n):
        temp = []
        for gene in human[i]:
            temp.append(genesClosest[allGenesList[random.randint(0,len(allGenesList))]])
        center = []
        for column in range(len(temp[0])):
            total = 0
            for vector in temp:
                total += vector[column]
            center.append(total/len(temp))
        cummulative = 0
        for entry in temp:
            cummulative += EuclidianDistance(entry,center)
        dist = cummulative/len(temp)
        if dist <= averageDist[i]:
            closer += 1
    pval.append(closer/n)
    print(str(closer/n))
    print(names[i] + " Finished!")

            
        

        
    
        


            
            







