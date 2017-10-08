def findOverlappingPair(str1, str2, str):
	max = -1
	len1 = len(str1)
	len2 = len(str2)

	for i in range(1, min(len1, len2)+1):
		if str1[len1-i:] == str2[:i]:
			if max < i:
				max = i
				str.clear()
				str.extend(str1 + str2[i:])

	for i in range(1, min(len1, len2)+1):
		if str1[:i] == str2[len2-i:]:
			if max < i:
				max = i
				str.clear()
				str.extend(str2 + str1[i:])

	return max

def findShortestSuperstring(arr):
	n_seqs = len(arr)
	while n_seqs != 1:
		max = -1
		l = 0
		r = 0
		result = []

		for i in range(0, n_seqs):
			for j in range(i+1, n_seqs):
				str = []
				res = findOverlappingPair(arr[i], arr[j], str)
				if max < res:
					max = res
					result = str
					l = i
					r = j
		n_seqs -= 1

		if max == -1:
			arr[0] += arr[n_seqs]
		else:
			arr[l] = result
			arr[r] = arr[n_seqs]


	return arr[0]

print(findShortestSuperstring([[1,2,3], [3,4], [1,2,3,4,5,6,7,8]]))
