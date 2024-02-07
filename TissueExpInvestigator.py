##Inestigate implication of codon usage on tissue expression
##Compare codon bias to tissue specific tRNA concentrations
##Use Tissue expression data to determine which genes are expressed in which tissues

##Processing input files
##Hash for converting Ensembl Identifies to UniProtKB
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

##Sorting Genes according to Tissue they are most highly expressed in

names = ["adipose", "adrenal", "marrow", "cortex", "colon", "duodenum", "endometrium", "esophagus", "fallopian", "gallbladder", "heart", "kidney", "liver",
         "lung", "lymph", "ovary", "pancreas", "placenta", "prostate", "rectum", "salivary", "skeletalmuscle", "smallintestine", "smoothmuscle", "spleen",
         "stomach", "testis", "thyroid", "tonsil", "bladder", "appendix", "skin"]
human = []

for i in names:
    human.append([])

for i in range(1,len(entries)):
    if entries[i][0] in key.keys() and key[entries[i][0]][0] != '':
        idtag = key[entries[i][0]][0]
        highest = 0.0
        tissue = 0
        for x in range(len(human)):
            if float(entries[i][x+2]) > highest:
                tissue = x
                highest = float(entries[i][x+2])
        human[tissue].append(idtag)


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
"""
##Calculate average distances
for tissue in human:
    scratch = []
    for gene in tissue:
        scratch.append(genes[gene])
         
                









    
