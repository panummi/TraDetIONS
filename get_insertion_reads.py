#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Script that takes a vcf and one sample at a time selects supporting reads for every insertion

import os
import sys, getopt
import pysam
import subprocess
from functools import lru_cache
from collections import defaultdict
from cyvcf2 import VCF
from pathlib import Path
from statistics import median


def count(v, sample_sam, window1, deviance, reads, samplesam, readnamefile, largest_read):	#supporting read selection
    insertion_size = v.format('LN')
    coordinates = str(v.format('CO')).split(',')[0]
    readcount_known= v.format('DR')
    readcount_known2 = 0
    for part in readcount_known:
            readcount_known2=part[1]+readcount_known2
    startco = coordinates.split('-')[0]
    endco = coordinates.split('-')[1].split("'")[0]
    sc = int(startco.split('_')[1])-window1
    ec = int(endco.split('_')[1])+window1
    window = insertion_size*deviance
    list_of_supplementary=[]
    supplementary_here=[]
    rows = pysam.view(sample_sam, '-f 2048  ', endco.split('_')[0] + ":"+ str(sc) + "-" + str(ec)).split('\n')
    insertion_sizes= {}
    for row in rows:
        name = row.split('\t')[0]
        if name in list_of_supplementary:
            supplementary_here.append(name)
        else:
            list_of_supplementary.append(name)
    for a in samplesam.fetch(endco.split('_')[0], sc, ec):
        bp_count = a.reference_start
        pos_in_read=1
        for ins in a.cigartuples:
            if ins[0] == 0:
                bp_count = bp_count + ins[1]
                pos_in_read=pos_in_read+ins[1]
            if ins[0] == 2:
                bp_count = bp_count + ins[1]
            if ins[0] == 1:
                if int(ins[1]) > insertion_size - window and int(ins[1]) < insertion_size + window and bp_count > sc and bp_count < ec:
                    insertion_sizes[a.qname] = [int(ins[1]), pos_in_read]
                    reads.write(a)
                    readnamefile.write(a.qname + "\n")
                    break
            if ins[0] == 4:
                pos_in_read=pos_in_read+ins[1]
            if ins[0] == 4 or ins[0] == 5:
                if bp_count >  int(startco.split('_')[1])-10 and bp_count < int(endco.split('_')[1])+10 and ins[1] > 30:
                    reads.write(a)
                    readnamefile.write(a.qname + "\n")
                    break
            if a.qname in supplementary_here:
                reads.write(a)
                readnamefile.write(a.qname + "\n")
                break

#selection of representative read

    sorted_d = dict(sorted(insertion_sizes.items(), key=lambda item: item[1]))
    length_of_dict = len(sorted_d)
    largest_read.write("type" + "\t" + "read_name" + "\t" +  "ins_size" + "\t" + "ins_pos_in_read" + "\n")

    if length_of_dict == 0:
        largest_read.write("largest" + "\t" + "." + "\t" +  "." + "\n")
        largest_read.write("median" + "\t" + "." + "\t" + "." + "\n")
    else:
        largest = (list(sorted_d.items())[length_of_dict-1])
        mediani = (list(sorted_d.items())[int((length_of_dict)/2)])
        largest_read.write("largest" + "\t" + largest[0] + "\t" +  str(largest[1][0]) + "\t" + str(largest[1][1]) + "\n")
        largest_read.write("median" + "\t" + mediani[0] + "\t" +  str(mediani[1][0]) + "\t" + str(mediani[1][1]) + "\n")


def main(argv):
    deviance = 0.5
    window = 200
    range="all"
    try:
        opts, args = getopt.getopt(argv, "hi:d:w:o:s:r:")
    except getopt.GetoptError:
        print("get_insertion_reads.py -h <help> -i <inputvcf> -d <allowed rate of deviance of length [0.5]> -w <window where possible reads are selected [200 bp]> -o <output directory> -s <sample directory> -r <range of samples to go through in the vcf [all]>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("get_insertion_reads.py -h <help> -i <inputvcf> -d <allowed rate of deviance of length [0.5]> -w <window where possible reads are selected [200 bp]> -o <output directory> -s <sample directory> -r <range of samples to go through in the vcf [all]>")
            sys.exit()
        elif opt in ("-d"):
            deviance = arg
        elif opt in ("-i"):
            inputfile = arg
        elif opt == ("-w"):
            window = arg
        elif opt == ("-o"):
            outputdir = arg
        elif opt == ("-s"):
            sampledir = arg
        elif opt == ("-r"):
            range = arg

    vcffile=VCF(inputfile)
    samples=vcffile.samples
    cont = False
    minimum = 0
    maximum = len(samples)-1
    if range != "all":
        minimum = int(range.split('-')[0])
        maximum = int(range.split('-')[1])
        print(minimum)
        print(maximum)
    u = 1
    for a in samples:	#go through the samples in vcf
        print(u)
        if u < minimum:
            u = u+1
            continue
        if u > maximum:
            break
        u = u+1
        print(a)
        name = (a[:-3])
#        pathfosample = name.split('/')[1].split('.')[0]
        Path(outputdir + "reads_in_insertions/" + a).mkdir(parents=True, exist_ok=True)
        sample_sam = ""
        with open(sampledir) as samples:
            for samp in samples:
                if a == samp.split("\t")[0]:
                    sample_sam=samp.split("\t")[1][:-1]
        if sample_sam == "":
            print("No cram file given for: " + a)
            break
        samplesam = pysam.AlignmentFile(sample_sam, "rb")
        reads = pysam.AlignmentFile(outputdir + a + ".bam", "wb", template=samplesam)
        vcf=VCF(inputfile, samples=a)

        for i, v in enumerate(vcf):	#go through insertions in sample
            if v.format('GT') != '' or v.format('AAL') != ['NAN']:
                vID = v.ID + "."  + "txt"
                largest_read = open(outputdir + "reads_in_insertions/" + a + "/" + vID + "_largest.txt", "w+")
                readnamefile = open(outputdir + "reads_in_insertions/" + a + "/" + vID, "w+")
                count(v, sample_sam, window, deviance, reads, samplesam, readnamefile, largest_read)

if __name__ == "__main__":
    main(sys.argv[1:])
