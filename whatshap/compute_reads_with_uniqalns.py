import sys
import random
import collections
from collections import defaultdict

# assumption ... all S' and before L's
filename = sys.argv[1]
f= open(filename, 'w')
#out = sys.argv[2]

d=defaultdict(set)

with open('all.x.sam') as fp:
	for line in fp:
		var=line.rstrip().split('\t')
		d[var[0]].add(var[1])
		
count=0
for k,v in d.items():
	if len(v) !=1:
		count= count+1
		f.write(str(k) + ',')

print(count)
