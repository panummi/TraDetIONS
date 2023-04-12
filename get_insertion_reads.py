#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys, getopt
import pysam
import subprocess
from functools import lru_cache
from collections import defaultdict
from cyvcf2 import VCF
from pathlib import Path
from statistics import median



def count(v, sample_sam, window1, deviance, reads, samplesam, readnamefile, largest_read):
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
                    print("HERE")
                    insertion_sizes[a.qname] = [int(ins[1]), pos_in_read] 
                    reads.write(a)
                    readnamefile.write(a.qname + "\n")
                    break
            if ins[0] == 4:
                pos_in_read=pos_in_read+ins[1]
            if ins[0] == 4 or ins[0] == 5:
                if bp_count >  int(startco.split('_')[1])-10 and bp_count < int(endco.split('_')[1])+10 and ins[1] > 30:
                    print("HERE")
                    reads.write(a)
                    readnamefile.write(a.qname + "\n")
                    break
            if a.qname in supplementary_here:
                print("HÄR")
                reads.write(a)
                readnamefile.write(a.qname + "\n")
                break
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
    
#


def main(argv):
    #cram = ""
    deviance = 0.5
    window = 200
    range="all"
    try:
        opts, args = getopt.getopt(argv, "hi:d:w:o:s:r:")
    except getopt.GetoptError:
        print("test.py [-S] -c <cram> -i <inputlist>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("test.py -c <cram> -i <inputputlist>")
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
    #readnamefile=open(outputdir + "/reads.txt", "w+")
    
    cont = False
    minimum = 0
    maximum = len(samples)-1
    if range != "all":
        minimum = int(range.split('-')[0])
        maximum = int(range.split('-')[1])
        print(minimum)
        print(maximum)
    u = 1
    for a in samples:
        print(u)
        if u < minimum:
            u = u+1
            continue    
        if u > maximum:
            break
        u = u+1
        print(a)
#        if a ==  "align/My_6112_T1_14_1706.sorted.map-ont.bam":
#            cont = True
#        if cont == False:
#            continue
        name = (a[:-3])
        pathfosample = name.split('/')[1].split('.')[0]
        #os.mkdir(outputdir + "reads_in_insertions/" + pathfosample) 
        Path(outputdir + "reads_in_insertions/" + pathfosample).mkdir(parents=True, exist_ok=True)
        namestart=a.split('/')[1].split('.')[0]
        sample_sam=sampledir+namestart+"/"+name+"cram"
        if a == "align/My_6304T1_16_0412.sorted.map-ont.bam":
            sample_sam = "/mnt/cgnano/projects/promethion/190827/ont_pipe/My_6304T1_16_0412_PAD92135/align/My_6304T1_16_0412.sorted.map-ont.cram"
        if a == "align/My_5006N_19_1461.sorted.map-ont.bam":
            sample_sam = "/mnt/cgnano/projects/promethion/ont_pipe_pass-current/My_5006N_19_1461_PAE46077/align/My_5006N_19_1461.sorted.map-ont.cram"
        if a == "align/My_5006T1_19_1460.sorted.map-ont.bam":
            sample_sam = "/mnt/cgnano/projects/promethion/ont_pipe_pass-current/My_5006T1_19_1460_PAE46514/align/My_5006T1_19_1460.sorted.map-ont.cram"
        if a == "align/My_6126T1_19_1307.sorted.map-ont.bam":
            sample_sam = "/mnt/cgnano/projects/promethion/ont_pipe_pass-current/My_6126T1_19_1307_PAE46577/align/My_6126T1_19_1307.sorted.map-ont.cram"
        if a == "align/My_8025T1_19_1458.sorted.map-ont.bam":
            sample_sam = "/mnt/cgnano/projects/promethion/ont_pipe_pass-current/My_8025T1_19_1458_PAE46517/align/My_8025T1_19_1458.sorted.map-ont.cram"
        if a == "align/My_6217N_19_1457.sorted.map-ont.bam":
            sample_sam = "/mnt/cgnano/projects/promethion/ont_pipe_pass-current/My_6217N_19_1457_PAE46448/align/My_6217N_19_1457.sorted.map-ont.cram"
        if a == "align/My_6143T1_19_1300.sorted.map-ont.bam":
            sample_sam = "/mnt/cgnano/projects/promethion/ont_pipe_pass-current/My_6143T1_19_1300_PAE46083/align/My_6143T1_19_1300.sorted.map-ont.cram"
        if a == "align/My_8025N_19_1459.sorted.map-ont.bam":
            sample_sam = "/mnt/cgnano/projects/promethion/ont_pipe_pass-current/My_8025N_19_1459_PAE46302/align/My_8025N_19_1459.sorted.map-ont.cram"
        if a == "align/c196_1_6833.sorted.map-ont.bam":
            sample_sam = "/mnt/cgnano/projects/promethion/ont_pipe_pass-current/c196_1_6833_PAE51989/align/c196_1_6833.sorted.map-ont.cram"
        if a == "align/c299_1_7246.sorted.map-ont.bam":
            sample_sam = "/mnt/cgnano/projects/promethion/ont_pipe_pass-current/c299_1_7246_PAE51990/align/c299_1_7246.sorted.map-ont.cram"
        if a == "align/c55_1_5391.sorted.map-ont.bam":
            sample_sam = "/mnt/cgnano/projects/promethion/ont_pipe_pass-current/c55_1_5391_PAE49467/align/c55_1_5391.sorted.map-ont.cram"
        if a == "align/c56_1_5395.sorted.map-ont.bam":
            sample_sam = "/mnt/cgnano/projects/promethion/ont_pipe_pass-current/c56_1_5395_PAE49640/align/c56_1_5395.sorted.map-ont.cram"
        if a == "align/s193_1_12_0763.sorted.map-ont.bam":
            sample_sam = "/mnt/cgnano/projects/promethion/ont_pipe_pass-current/s193_1_12_0763_PAE49478/align/s193_1_12_0763.sorted.map-ont.cram"
        if a == "align/My_6091T1_14_1017.sorted.map-ont.bam":
            sample_sam = "/mnt/cgnano/projects/promethion/ont_pipe_pass-current/My_6091T1_14_1017_PAD64263/align/My_6091T1_14_1017.sorted.map-ont.cram"
        if a == "align/Fam_c596_1_8259TK.sorted.map-ont.bam":
            sample_sam = "/mnt/cgnano/projects/promethion/ont_pipe_pass-current/Fam_c596_1_8259TK_PAD65691/align/Fam_c596_1_8259TK.sorted.map-ont.cram"
        samplesam = pysam.AlignmentFile(sample_sam, "rb")
        reads = pysam.AlignmentFile(outputdir + name.split('/')[1] + "bam", "wb", template=samplesam)
        vcf=VCF(inputfile, samples=a)
        print(a)
#        print(vcf)
        for i, v in enumerate(vcf):
            if v.format('GT') != '' or v.format('AAL') != ['NAN']:
                vID = v.ID + "."  + "txt"
                #if vID == "S_1.txt":
                print(vID)
                largest_read = open(outputdir + "reads_in_insertions/" + pathfosample+ "/" + vID + "_largest.txt", "w+")
                readnamefile = open(outputdir + "reads_in_insertions/" + pathfosample+ "/" + vID, "w+")
#                readnamefile.write(vID + "\n")
                count(v, sample_sam, window, deviance, reads, samplesam, readnamefile, largest_read)


if __name__ == "__main__":
    main(sys.argv[1:])
