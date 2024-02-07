import pandas as pd


ref2Ens = {}

with open("Ensembl-RefSeq-UniProt.csv","r") as file:
    data = file.readlines()
    for i in data:
        line = i.strip().split(',')
        if line[5] != '':
            ref2Ens[line[5]] = line[0]

protConc = {}

with open("TissueProteinLevels.csv", "r") as file:
    data = file.readlines()
    names = data[0].strip().split(',')[1:]
    for line in data[1:]:
        temp = line.strip().split(',')
        iD = temp[0].split('.')[0]
        if iD in ref2Ens.keys():
            protConc[ref2Ens[iD]] = temp[1:]

geneExp = {}

with open("TissueExp.csv", "r") as file:
    data = file.readlines()
    for line in data[1:]:
        temp = line.strip().split(',')
        geneExp[temp[0]] = temp[1:]

codonBias = {}

with open("Human101.UNIQ.codon_bias.csv") as file:
    data = file.readlines()
    codons = data[0].strip().split(',')[1:]
    for line in data[1:]:
        temp = line.strip().split(',')
        codonBias[temp[0].split('|')[0]] = temp[1:]

#intersection of 3 sets

pc = set(protConc.keys())
ge = set(geneExp.keys())
cb = set(codonBias.keys())

intersection = list(pc & ge & cb)

for key in pc:
    if key not in intersection:
        del protConc[key]

for key in ge:
    if key not in intersection:
        del geneExp[key]

for key in cb:
    if key not in intersection:
        del codonBias[key]

normLvls = {}
for gene in geneExp.keys():
    for tissue in range(len(names)):
        protConc[gene][tissue] = float(protConc[gene][tissue])
        geneExp[gene][tissue] = float(geneExp[gene][tissue])
    expTotal = sum(geneExp[gene])
    protTotal = sum(protConc[gene])
    rLvls = []
    for tissue in range(len(names)):
        exp = float(geneExp[gene][tissue])
        prot = protConc[gene][tissue]
        if expTotal == 0:
            rLvls.append(str(0))
        elif protTotal == 0:
            rLvls.append(str(0))
        elif exp == 0:
            rLvls.append(str(0))
        else:
            level = (prot/protTotal)/(exp/expTotal)
            if level <= 1.5:
                rLvls.append(str(0))
            else:
                rLvls.append(str(1))

    normLvls[gene] = rLvls

with open("TissueExpConcRatio1.5.csv", "w") as file:
    print("Gene ID," + ",".join(names),file = file)
    for gene in normLvls.keys():
        print(gene + "," + ",".join(normLvls[gene]), file = file)

with open("Human101TisGenes.csv", "w") as file:
    print("Gene ID," + ",".join(codons), file = file)
    for gene in normLvls.keys():
        print(gene + "," + ",".join(codonBias[gene]), file = file)

with open("TisTRatio1.5Classes.csv", "w") as file:
    for gene in normLvls.keys():
        count = 0
        for i in range(len(normLvls[gene])):
            if int(normLvls[gene][i]) == 1:
                print(str(i) + "," + ",".join(codonBias[gene]), file = file)
                count += 1
        if count == 0:
            print(str(len(names)) + "," + ",".join(codonBias[gene]), file = file)

df = pd.read_csv("TisTRatio1.5Classes.csv", header = None)

df = df.drop_duplicates()

df = df.sort_values(df.columns[0])

df.to_csv("TisTRatio1.5Classes.csv", header = False, index = False)
        
    
        
    
    

    
    
    
        
    
    



