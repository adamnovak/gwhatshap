import sys

# assumption ... all S' and before L's
filename = sys.argv[1]
#out = sys.argv[2]

d={}
count=1

with open(filename) as fp:
	for line in fp:
		var=line.split(	)
		if(var[0] == 'L'):
			if var[2] == '+' and var[4].rstrip() == '+':
				d[var[1]]=count
				count+=1
				d[var[3]]=count
				count+=1
			if var[2] == '-' and var[4].rstrip() == '+':
				tmp1 = 2*int(var[1])
				d[str(tmp1)]=count
				count+=1
				d[var[3]]=count
				count+=1
			if var[2] == '+' and var[4].rstrip() == '-':
				tmp1 = 2*int(var[3])
				d[var[1]]=count
				count+=1
				d[str(tmp1)]=count
				count+=1
			if var[2] == '-' and var[4].rstrip() == '-':
				tmp1 =  2*int(var[1])
				tmp2 =  2*int(var[3])
				d[str(tmp1)] = count
				count = count + 1
				d[str(tmp2)] = count
				count = count + 1

with open(filename) as fp:
	for line in fp:
		var=line.split(	)
		if(var[0] == 'L'):
			if var[2] == '+' and var[4].rstrip() == '+':
				print(str(d[var[1]])+ " " + str(d[var[3]]))
			if var[2] == '-' and var[4].rstrip() == '+':
				tmp1 = 2*int(var[1])
				print(str(d[str(tmp1)])+ " " + str(d[var[3]]))
			if var[2] == '+' and var[4].rstrip() == '-':
				tmp1 = 2*int(var[3])
				print(str(d[var[1]])+ " " + str(d[str(tmp1)]))
			if var[2] == '-' and var[4].rstrip() == '-':
				tmp1 =  2*int(var[1])
				tmp2 =  2*int(var[3])
				print(str(d[str(tmp1)])+ " " + str(d[str(tmp2)]))

print(count)

