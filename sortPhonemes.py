from operator import itemgetter

d=[]
with open("labelFile06BASIC.lab", "r") as infile:
    for line in infile.readlines():
        data = line.strip().split(' ')
        data = [float(data[0]) , float(data[1]), (data[2])]
        d.append(data) # add the line to a list `d`
#for line in sorted(lines, key=lambda line: line.split()[0]):
d.sort(key = itemgetter(0))
print(d)

with open("sortedLabelBasicWihSilence06", 'w+') as outfile:
    for line in d:
        outfile.write(str(line[0]) + " " + str(line[1]) + " "+ line[2]+"\n");
        #outfile.write(''.join(str(line))+'\n');
