import vcf
import sys

input_file = sys.argv[1]

vcf_reader = vcf.Reader(open(input_file, 'r', encoding='utf-8'))
block_ids = set()

for record in vcf_reader:
    try:
        block_ids.add(record.samples[0]['PS'])
    except:
        continue

#print(len(block_ids))

for i in block_ids:
    vcf_reader = vcf.Reader(open(input_file, 'r', encoding='utf-8'))
    vcf_writer = vcf.Writer(open(input_file+ str(i) + '.vcf', 'w'), vcf_reader)
    for record in vcf_reader:
        try:
            #print(record.samples[0]['PS'])
            #print(str(i))
            if str(record.samples[0]['PS']) == str(i):
                vcf_writer.write_record(record)
        except:
            continue 





