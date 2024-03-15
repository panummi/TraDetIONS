
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



def main(argv):
    #cram = ""
    deviance = 0.5
    window = 200
    range="all"
    try:
        opts, args = getopt.getopt(argv, "hi:d:w:o:s:r:")
    except getopt.GetoptError:
        print("potential_transductions.py -h <help> -i <inputvcf> -d <allowed rate of deviance of length [0.5]> -w <window where possible reads are selected [200 bp]> -o <output directory> -s <sample directory> -r <range of samples to go through in the vcf [all]>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("potential_transductions.py -h <help> -i <inputvcf> -d <allowed rate of deviance of length [0.5]> -w <window where possible reads are selected [200 bp]> -o <output directory> -s <sample directory> -r <range of samples to go through in the vcf [all]>")
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

        print(a)
#        print(vcf)
        for i, v in enumerate(vcf):
            if v.format('GT') != '':
                vID = v.ID + "."  + "txt"
                #if vID == "S_1.txt":
                print(vID)
                largest_read = open(outputdir + "reads_in_insertions/" + a + "/" + vID + "_largest.txt", "w+")
                readnamefile = open(outputdir + "reads_in_insertions/" + a + "/" + vID, "w+")
#                readnamefile.write(vID + "\n")
                count(v, sample_sam, window, deviance, reads, samplesam, readnamefile, largest_read)


if __name__ == "__main__":
    main(sys.argv[1:])

