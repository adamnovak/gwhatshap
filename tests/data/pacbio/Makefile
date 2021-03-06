# This Makefile documents how the PacBio test data set was generated.
# The data is from the Genome in a Bottle dataset (individual HG004)

SHELL := /bin/bash

all: reference.fasta pacbio.bam pacbio.bam.bai variants.vcf hapcut.txt

ucsc.hg19.fasta:
	@echo "Please obtain file ucsc.hg19.fasta from the GATK resource bundle to continue"
	@exit 1

reference.fasta: ucsc.hg19.fasta ucsc.hg19.fasta.fai
	( echo '>ref'; samtools faidx $< chr6:10029001-10055081 | sed 1d | tr a-z A-Z ) > $@

downloaded-illumina.bam:
	samtools view -b ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/NIST_Illumina_2x250bps/HG004_250bp_All-29933060/HG004-250bp-All_S1.bam chr6:10039500-10049000 > $@

downloaded-pacbio.bam:
	samtools view -b ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/PacBio_MtSinai_NIST/CSHL_bwamem_bam_GRCh37/BWA-MEM_Chr6_HG004_merged_11_12.sort.bam 6:10039500-10049000 > $@

%.bam: downloaded-%.bam
	( echo -e "@HD\tVN:1.2\tSO:coordinate\n@SQ\tSN:ref\tLN:27000\n@RG\tID:1\tSM:HG004_250bp_All"; samtools view -s 0.3 -x MD -x RG $< | awk -vFS='\t' -vOFS='\t' '{$$3="ref"; $$4-=10029000; print $$0 "\tRG:Z:1"}' ) | samtools view -b - > $@

illumina-full.bam: downloaded-illumina.bam
	( echo -e "@HD\tVN:1.2\tSO:coordinate\n@SQ\tSN:ref\tLN:27000\n@RG\tID:1\tSM:HG004_250bp_All"; samtools view -x MD -x RG $< | awk -vFS='\t' -vOFS='\t' '{$$3="ref"; $$4-=10029000; print $$0 "\tRG:Z:1"}' ) | samtools view -b - > $@

illumina.bam: illumina-full.bam
	samtools view -b -s 0.3 $< > $@

variants.vcf: reference.fasta illumina-full.bam illumina.bam.bai
	freebayes -f reference.fasta illumina-full.bam | \
		awk -vFS='\t' -vOFS='\t' '!/^#/ {$$8="."};1' | \
		grep -v '^##INFO=' | \
		sed -r 's|GT(:[A-Z][A-Z])*\t([01/]*):.*|GT\t\2|' > $@ && \
	echo -e 'ref\t26081\t.\tG\tT\t17.3\t.\t.\tGT\t0/1' >> $@

hairs.txt: reference.fasta variants.vcf pacbio.bam pacbio.bam.bai
	extractHAIRS --ref reference.fasta --indels 1 --VCF variants.vcf --bam pacbio.bam --maxIS 600 > $@

hapcut.txt: hairs.txt variants.vcf
	HAPCUT --fragments hairs.txt --VCF variants.vcf --output hapcut.txt

%.bam.bai: %.bam
	samtools index $<

%.fasta.fai: %.fasta
	samtools faidx $<
