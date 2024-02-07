## Testing point packing for codon selection within a single amino acid
## Explore the relationship between the resolution of the point space and
# the minimum distance of the point packer.

#Resolutoin Dial
rez = 3

plane_three = []
for x in range(rez + 1):
    for y in range((rez-x) + 1):
        z = rez - (x + y)
        plane_three.append([x/rez, y/rez, z/rez])

#Euclidian distance metric
def EDist(point,cluster):
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

#Evolutionary Point Packing Algorithm
#libraries
from numpy import random
from random import shuffle
simplex = plane_three
mind = float(47/100)
ecosystem = []
rmr = 30
size = 40
numofmatings = 30
generations = 30
goats = []

for i in range(size):
    container = []
    for i in range(rmr):
        sample = random.randint(0,len(simplex))
        container.append(simplex[sample])
    pointpack = []
    for i in container:
        if len(pointpack) == 0:
            pointpack.append(i)
        else:
            x = 0
            tag = True
            while x < len(pointpack) and tag == True:
                if EDist(i,pointpack[x]) > mind:
                    x += 1
                else:
                    tag = False
            if tag == True:
                pointpack.append(i)
    ecosystem.append(pointpack)


for i in range(generations):
    ecosystem.sort(key = lambda s: len(s), reverse = True)
    for i in ecosystem[:3]:
        goats.append(i)
    nextgen = ecosystem
    for i in range(numofmatings):
        female = random.randint(0,len(ecosystem))
        male = random.randint(0,len(ecosystem))
        while male == female:
            male = random.randint(0,len(ecosystem))
        genepool = []
        for i in ecosystem[male]:
            genepool.append(i)
        for i in ecosystem[female]:
            genepool.append(i)
        for i in range(rmr):
            genepool.append(simplex[random.randint(0,len(simplex))])
        shuffle(genepool)
        child = []
        for i in genepool:
            if len(child) == 0:
                child.append(i)
            else:
                x = 0
                tag = True
                while x < len(child) and tag == True:
                    if EDist(i,child[x]) > mind:
                        x += 1
                    else:
                        tag = False
                if tag == True:
                    child.append(i)
        spot = random.randint(0,len(nextgen))
        if len(child) >= len(nextgen[spot]):
            nextgen[spot] = child
    ecosystem = nextgen




            
