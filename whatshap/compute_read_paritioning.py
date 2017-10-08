from collections import defaultdict

def compute_read_partitioning_accuracy():
	true_hap1 =[]
	true_hap2= []
	#f = open(true_file, 'r')
	
	for line in open('true_partioning'):
		var = line.rstrip().split(" ")
		#print(var)
		if var[1] =='1':
			true_hap1.append(var[0])
		else:
			true_hap2.append(var[0])
	#print(len(true_hap1))
	#print(true_hap1)
		
	pred_hap1 = defaultdict(list)
	pred_hap2 = defaultdict(list)
	total=0
	blocks =set()
	for line in open('whole_genome.predicted_read_partionting.pred_old'):
		tokens= line.rstrip().split(" ")
		blocks.add(tokens[1])
		if tokens[2] =='1':
			pred_hap1[tokens[1]].append(tokens[0])
		else:
			pred_hap2[tokens[1]].append(tokens[0])
		total+=1
	print('predicted')
	#print(pred_hap1['6878'])
	#print(pred_hap2['6878'])
	
	print(len(list(blocks)))
	count = 0
	#duplicates = [24899,18733,3763,5241,18841,17801,4623,17491,12705,1679,13699,18466,19418,24090,23744,13931,15366,13388,6334,2,8913,10192,20615,24999,24254,9050,15084,24061,5556,24708,10257,23275,18922,19628,23902,8623,22461,13405,21608,4279,17950,15274,15598]
	for k in blocks:
	#	if k in duplicates:
	#		continue
		count = count+max(len(list(set(pred_hap1[k]).intersection(set(true_hap1)))), len(list(set(pred_hap1[k]).intersection(set(true_hap2)))))
		#print(k, count)
		count= count+ max(len(list(set(pred_hap2[k]).intersection(set(true_hap1)))), len(list(set(pred_hap2[k]).intersection(set(true_hap2)))))
		#print(k, count)
		#print(k)
	percent_partitionining_accuracy =  count/total
	print('percent_partitionining_accuracy')
	print(percent_partitionining_accuracy)

compute_read_partitioning_accuracy()
