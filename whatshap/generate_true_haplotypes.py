import sys
from collections import defaultdict

variants_file = sys.argv[1]
nodes_in_contigs_file = sys.argv[2] 
out_file=sys.argv[3]

output= open(out_file, 'w')
locus_branch_mapping_tmp=defaultdict(list)

with open(variants_file) as fp:
	for line in fp:
		if line[0]==">":
			var=line.split("/")
			bubbleID = int(var[0][9:])
			tmpvector = line.split(" ")
			branchID = int(var[1].split(" ")[0])
			tmp=list()
			for n in range(len(tmpvector)-2,0,-1):
				if tmpvector[n]!="VertexIDs:":
					tmp.append(tmpvector[n])
				else:
					break
			locus_branch_mapping_tmp[bubbleID].append(tmp)

#print(locus_branch_mapping_tmp)


bubble_mapping=defaultdict()
for key,value in locus_branch_mapping_tmp.items():
	if len(value) ==2:
		if len(list(value[0]))==1 and len(list(value[1]))==1:
			str1=" ".join(str(x) for x in value[0])
			str2 = " ".join(str(x) for x in value[1])
			bubble_mapping[str1] = str2

#print(bubble_mapping)


with open(nodes_in_contigs_file) as fp:
	for line in fp:
		print(line)
		strprint = ""
		printtmp=[]
		line = line.lstrip()
		var= line.split(" ")
		#print(var)
		for i in var:
			if i in bubble_mapping:
				printtmp.append(bubble_mapping[i])
			else:
				printtmp.append(i)
		strprint = " ".join(x for x in printtmp)
		output.write(" "+ strprint)
		



