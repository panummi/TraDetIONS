#!/usr/bin/env python
# -*- coding: utf-8 -*-

#goes thorugh vcf and marks the number of samples sharing the insertions and recognizes if the insertion is somatic

import os
import sys, getopt
import pysam


def main(argv):
	vcf_file=''
	output_file=''
	sample_file=''
	try:
		opts, args = getopt.getopt(argv,"hv:o:s:")
	except getopt.GetoptError:
		print ('mark_variants.py -v <vcf-input> -o <outputfile> -s <sample>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print ('mark_variants.py -v <vcf-input> -o <outputfile> -s <sample>')
			sys.exit()
		elif opt in ("-v"):
			vcf_file = arg
		elif opt in ("-o"):
			output_file = arg
		elif opt in ("-s"):
			sample_file = arg

	def is_somatic(variant, samples):
		vector=[s['GT'] for s in variant.samples.values()]
		NA=0
		NG=0
		NU=0
		S=0
		TA=0
		TG=0
		for a, ch in enumerate(vector, 0):
			matching=""
			sample=samples[a].split('\t')[0]
			if ch != (None, None):
				pairid=samples[a].split('\t')[2]
				if pairid != '-':
					if vector[int(pairid)] != (None, None):
						matching="match"
					else:
						matching="no"
				else:
					matching="absent"
				if samples[a].split('\t')[1] == "N":
					if matching=="no":
						variant.samples[sample]['STATUS']="Normal_unique"
						NU += 1
					if matching=="match":
						variant.samples[sample]['STATUS']="Normal_germline"
						NG += 1
					if matching=="absent":
						variant.samples[sample]['STATUS']="Normal_absent_pair"
						NA += 1
				if samples[a].split('\t')[1] == "T":
					if matching=="no":
						variant.samples[sample]['STATUS']="Somatic"
						S += 1
					if matching=="match":
						variant.samples[sample]['STATUS']="Tumor_germline"
						TG += 1
					if matching=="absent":
						variant.samples[sample]['STATUS']="Tumor_absent_pair"
						TA+=1
			else:
				variant.samples[sample]['STATUS']="NA"
			variant.info['NU']=NU
			variant.info['NG']=NG
			variant.info['NA']=NA
			variant.info['S']=S
			variant.info['TG']=TG
			variant.info['TA']=TA

	vcf = pysam.VariantFile(vcf_file)

	samples = open(sample_file).readlines()

	vcf.header.formats.add("STATUS",".","String","Is variant somatic/germline")
	vcf.header.info.add('NU', number=1, type='Integer', description='Number of times the varinat is normal_unique in samples')
	vcf.header.info.add('NG', number=1, type='Integer', description='Number of times the variant is normal_germline in samples')
	vcf.header.info.add('NA', number=1, type='Integer', description='Number of times the variant is normal_absent_pair in samples')
	vcf.header.info.add('S', number=1, type='Integer', description='Number of times the variant is somatic in samples')
	vcf.header.info.add('TG', number=1, type='Integer', description='Number of times the variant is tumor_germline in samples')
	vcf.header.info.add('TA', number=1, type='Integer', description='Number of times the variant is tumor_absent_pair in samples')
	output=pysam.VariantFile(output_file, "w", header=vcf.header)

	for v in vcf:
		is_somatic(v, samples)
		output.write(v)
	vcf.close()

if __name__ == "__main__":
	main(sys.argv[1:])
