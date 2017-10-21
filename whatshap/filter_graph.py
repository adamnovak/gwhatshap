import sys

# assumption ... all S' and before L's
filename = sys.argv[1]
#out = sys.argv[2]

ids = [978511, 978512, 978513, 978514, 978515, 978516, 978517, 978518, 978519, 978520, 978521, 978522, 978523, 978524, 978525, 978526, 978527, 978528, 978529, 978530, 978531, 978532, 978533, 978534, 978535, 978536, 978537, 978538, 978539, 978540, 978541, 978542, 978543, 978544, 978545, 978546, 978547, 978550, 978551, 978552, 978553, 978554, 978555, 978556, 978557, 978558, 978576, 978582, 978707, 978708, 978709, 978710, 978711, 978712, 978713, 978714, 978715, 978716, 978717, 978954, 979139, 979140, 979460, 979879, 979880, 979995, 979996, 979997, 981065, 981526, 981927, 982844, 982845, 982846, 982847, 982848, 983271, 983960, 984155, 984157, 984158, 984159, 985456, 985698, 986649, 988125, 988391, 988839, 989173, 989459, 989460, 989461, 989508, 990563, 990564, 990982, 991334, 991335, 991336, 991337, 992250, 992847, 992848, 993248, 993463, 993855, 994559, 994688]


with open(filename) as fp:
	for line in fp:
		var=line.split('\t')
		if(var[0] == 'S') and int(var[1]) not in ids:
			print(line)
		if(var[0] == 'L') and int(var[1]) not in ids and (var[3]) not in ids:
		else:
			print(line)
