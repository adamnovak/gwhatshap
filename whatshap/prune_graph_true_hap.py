import sys
from networkx import nx

nodes_in_contigs_file=sys.argv[1]
filename = sys.argv[2]

out = sys.argv[3]

output= open(out, 'w')

tmp1=[]
tmp2=[]
with open(nodes_in_contigs_file) as fp:
	for line in fp:
		line = line.lstrip()
		var= line.split(" ")
		for i in range(0,len(var),2):
			tmp1.append(var[i])
		for i in range(1,len(var),2):
			tmp2.append(var[i])


with open(filename) as fp:
	for line in fp:
		var=line.split("\t")
		if(var[0] == 'S'):
			if var[1] in tmp1 or var[1] in tmp2:
				output.write(line)
		else:
			if var[1] in tmp2 and var[3] in tmp1:
				output.write(line)
			





