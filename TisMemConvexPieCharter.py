import matplotlib.pyplot as plt
import os

script_dir = os.path.dirname(__file__)
results_dir = os.path.join(script_dir, 'Tissue Membership Pie Charts/')


folder = ["Enriched", "Enhanced", "Elevated"]

protocol = folder[2]

index = {}

with open(protocol + "TisMem.csv","r") as file:
    data = file.readlines()
    treatment = data[0].strip().split(',')[1:]
    for line in data[1:]:
        temp = line.strip().split(',')
        if temp[0] in index.keys():
            index[temp[0]].append(temp[1:])
        else:
            index[temp[0]] = [temp[1:]]

background = []
labels = []
for tissue in index.keys():
    count = 0
    for gene in index[tissue]:
        count +=1
    background.append(count)
    labels.append(tissue)



for i in range(len(treatment)):
    mem = []
    for tissue in index.keys():
        count = 0
        for gene in index[tissue]:
            if gene[i] == "1":
                count += 1
        mem.append(count)
    fig, axs = plt.subplots(2)
    plt.figure(figsize = (12.0,9.0))
    axs[0].pie(background, labels=labels, labeldistance = None, autopct='%1.1f%%', pctdistance = 1.2, textprops = {'fontsize': 2.5})
    axs[0].set_title("Backgound")
    axs[1].pie(mem, labels=labels, labeldistance = None, autopct='%1.1f%%', pctdistance = 1.3, textprops = {'fontsize': 2.5})
    axs[1].set_title("Convex Hull")
    fig.suptitle(protocol + " Expression " + labels[i//3] + " Membership Breakdown")
    fig.legend(labels, loc = "center left", prop = {"size":6})
    fig.savefig(fname = results_dir + protocol + " Convex/" + str(labels[i//3]), dpi = 900, bbox_inches = "tight")
    fig.clf()
    
    
        

"""    
plt.figure(figsize=(9.0,7.0))
plt.scatter(x1, y1, c = z1, s = 1, cmap=plt.cm.RdYlGn, vmax = upper, vmin = lower)
plt.xlabel("\nX (m)")
plt.ylabel("Y (m)\n")
plt.title("Profit Mapping for " + farmname +" Farm\n" + crop + "  " + year + "\n")
plt.colorbar().ax.set_ylabel("Profitability ($/Acre)")
plt.savefig(fname = results_dir + farmname + "_" + crop + "_" + year, dpi = 900, bbox_inches = "tight", pad_inches = 0.5)
"""
