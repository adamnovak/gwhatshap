import vcf
import sys

# get haplotigs based on vcf and canu contigs
input_file = sys.argv[1]
canu_contig_file = sys.argv[2]

vcf_reader = vcf.Reader(open(input_file, 'r', encoding='utf-8'))
block_ids = set()

for record in vcf_reader:
    try:
        block_ids.add(record.samples[0]['PS'])
    except:
        continue

canu_contig_ref = list()
sorted_block_ids = sorted(block_ids)
with open(canu_contig_file) as fp:
	for line in fp:
		var=line.rstrip()[0]
		if var!='>':
                    for x in list(line.rstrip()):
                        canu_contig_ref.append(x)

canu_contig_alt = list()
for a in canu_contig_ref:
	canu_contig_alt.append(a)

#print(len(block_ids))
#print(canu_contig_alt)
for i in range(0,len(sorted_block_ids)):
    vcf_reader = vcf.Reader(open(input_file, 'r', encoding='utf-8'))
    vcf_writer = open(input_file+ str(sorted_block_ids[i]) + '.fasta', 'w')
    #vcf_writer2 = open(input_file+ str(sorted_block_ids[i]) + '.fasta', 'w')
    hap1_seq = ''
    hap2_seq = ''
    if i==0:
        canu_contig_alt_tmp = canu_contig_alt[0:sorted_block_ids[i]-1]
        canu_contig_ref_tmp = canu_contig_ref[0:sorted_block_ids[i]-1]
    else:
        canu_contig_alt_tmp = canu_contig_alt[sorted_block_ids[i-1]-1:sorted_block_ids[i]-1]
        canu_contig_ref_tmp = canu_contig_ref[sorted_block_ids[i-1]-1:sorted_block_ids[i]-1]
    #print(canu_contig_alt_tmp)
    for record in vcf_reader:
        try:
            #print(record.samples[0]['PS'])
            #print(str(i))

            if str(record.samples[0]['PS']) == str(sorted_block_ids[i]):
                pos=int(record.POS)
                ref=record.REF
                alt=record.ALT[0]
                #print(alt)
                hap1=record.samples[0]['GT']
                allele1=hap1.split('|')
                if int(allele1[0])==0:
                    canu_contig_ref_tmp[pos]=str(ref)
                else:
                    canu_contig_ref_tmp[pos]=str(alt)
    
                if int(allele1[1])==0:
                    canu_contig_alt_tmp[pos]=str(ref)
                else:
                    canu_contig_alt_tmp[pos]=str(alt)                
        except:
            continue 
        #print(canu_contig_ref_tmp)
    for p in  canu_contig_ref_tmp:
        hap1_seq = hap1_seq+p
    for q in  canu_contig_alt_tmp:
        hap2_seq = hap2_seq+q
    vcf_writer.write(">" + input_file+ str(sorted_block_ids[i]) + "_1" + "\n")
    vcf_writer.write(hap1_seq + "\n")
    vcf_writer.write(">" + input_file+ str(sorted_block_ids[i]) + "_2" + "\n")
    vcf_writer.write(hap2_seq + "\n")






