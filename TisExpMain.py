def EuclidianDistance(point,cluster):
    if len(point) != len(cluster):
        return False
    else:
        tot = 0
        for i in range(len(point)):
            tot += ((float(point[i])-float(cluster[i]))**2)                           
        dist = (tot)**0.5
        return dist

def main():
    genes = {}
    with open("Homo_sapiens.GRCh38.codon_bias.csv", "r") as file:
        data = file.readlines()
        for i in data[1:]:
            line = i.strip().split(',')
            name = line[0]
            bias = line[1:]
            genes[name] = bias

    data = []
    for i in genes.keys():
        temp = []
        for col in genes[i]:
            temp.append(float(col))
        data.append(temp)

    print("Assembled")

    total = 0
    n = 0
    for i in range(len(data)):
        for j in range(i+1,len(data)):
            total += EuclidianDistance(data[i],data[j])
            n += 1

    mu = total/n

    sqtotal = 0
    for i in range(len(data)):
        for j in range(i+1,len(data)):
            sqtotal += ((EuclidianDistance(data[i],data[j]) - mu)**2)

    var = sqtotal/n
    sd = var**0.5

main()
