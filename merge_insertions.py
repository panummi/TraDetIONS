#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Script merges insertions in a vcf file

import os
import sys, getopt
from typing import Sequence
import pysam
import subprocess
from functools import lru_cache
from collections import defaultdict
from pathlib import Path
from collections import namedtuple
import mappy
import logging as log
from operator import attrgetter
from operator import itemgetter


def get_info(v, TE_single, TE_cluster):	#returns info about individual insertion and insertion cluster
    single = TE_single(v.chrom, v.pos, len(v.alts[0]), v.id, v)
    cluster = TE_cluster(v.chrom, v.pos, v.pos, len(v.alts[0]), len(v.alts[0]), [v.id], [v], [single])
    return single, cluster

def matchnumbers(a, b1, b2, limit, comparison):	#check if insertions fill all merging criteria
    if (a > b1 and a < b2):
        return True
    if comparison == "pos" and (abs(a - b2) < limit) or (abs(a - b1) < limit):
        return True
    if comparison == "length" and (abs((a-(b1))/abs(b2)) < limit or abs((a-(b1))/abs(a)) < limit):
        return True
    else:
        return False

def comparenear(a,b):	#check if insertions are close
    if a.chr == b.chr and matchnumbers(a.pos,b.minpos,b.maxpos, 40, "pos"):
        return True
    else:
        return False

def compare(a, b):	#check if insertions lengths match
    if (a.length,b.minlength,b.maxlength, 0.75, "length") and a.variant_id not in b.variant_ids:
        return True
    else:
        return False

def make_cluster(a, cluster, TE_cluster):	#makes cluster of insertions
    te_id_list = cluster.variant_ids.copy()
    te_id_list.append(a.variant_id)
    te_list = cluster.variants.copy()
    te_list.append(a.variant)
    single_list=cluster.TE_singles.copy()
    single_list.append(a)
    if a.pos < cluster.minpos:
        minpos = a.pos
    else:
        minpos = cluster.minpos
    if a.pos > cluster.maxpos:
        maxpos = a.pos
    else:
        maxpos = cluster.maxpos
    if a.length < cluster.minlength:
        minlength = a.length
    else:
        minlength = cluster.minlength
    if a.length > cluster.maxlength:
        maxlength = a.length
    else:
        maxlength = cluster.maxlength
    return(TE_cluster(a.chr, minpos, maxpos, minlength, maxlength, te_id_list, te_list, single_list))

def merge_output(cluster, germl):	#create output for merged variant cluster
    id="P_"+str(germl)
    germl=germl+1
    sum=0
    i=0
    TES=''
    bestmapq=-1
    best_seq=''
    variant=cluster.TE_singles[0].variant.copy()
    NU=0
    NG=0
    NA=0
    S=0
    TG=0
    TA=0

    variant.info.__delitem__("SUPP")
    variant.info.__delitem__("SUPP_VEC")
    variant.info.__delitem__("SVLEN")
    variant.info.__delitem__("SVTYPE")
    variant.info.__delitem__("SVMETHOD")
    variant.info.__delitem__("CIPOS")
    variant.info.__delitem__("CIEND")
    variant.info.__delitem__("STRANDS")
    if "POSITION" in variant.info:
        variant.info.__delitem__("POSITION")
    if "FLANKINGstart" in variant.info:
        variant.info.__delitem__("FLANKINGstart")
    if "FLANKINGend" in variant.info:
        variant.info.__delitem__("FLANKINGend")
    if "SVclass" in variant.info:
        variant.info.__delitem__("SVclass")
    if "poly" in variant.info:
        variant.info.__delitem__("poly")

    for a in cluster.TE_singles:
        v=a.variant
        sum=v.pos+sum
        i=i+1
        endsum=0
        TES=TES+v.id+":"
        best_seq=v.alts[0]
        NU=v.info['NU']+NU
        NG=v.info['NG']+NG
        NA=v.info['NA']+NA
        S=v.info['S']+S
        TG=v.info['TG']+TG
        TA=v.info['TA']+TA
        for s in v.samples:
            if v.samples[s] != variant.samples[s]:
                if v.samples[s]['GT'] != (None, None):
                    for w in v.samples[s]:
                        if type(variant.samples[s][w]) == str and  type(v.samples[s][w]) == tuple:
                            variant.samples[s][w] = v.samples[s][w][0]
                        else:
                            variant.samples[s][w] = v.samples[s][w]
        if v.stop:
            endsum=v.stop + endsum

    variant.id = id
    variant.pos=int(sum/i)
    variant.alleles=(['N', best_seq])
    variant.stop = int(endsum/i)
    variant.info['NU'] = NU
    variant.info['NG'] = NG
    variant.info['NA'] = NA
    variant.info['S'] = S
    variant.info['TG'] = TG
    variant.info['TA'] = TA

    return(variant, germl, [id, TES[:-1]])


def main(argv):
    outputfile=''
    inputl=''
    try:
        opts, args = getopt.getopt(argv, "hi:o:")
    except getopt.GetoptError:
        print("merge_insertion.py -i <inputvcf> -o <outputvcf>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("merge_insertion.py -i <inputvcf> -o <outputvcf>")
            sys.exit()
        elif opt in ("-i"):
            inputl = arg
        elif opt in ("-o"):
            outputfile = arg

    list_of_singles=[]
    list_of_clusters=[]
    TE_single = namedtuple('TE_single', 'chr pos length variant_id variant')
    TE_cluster = namedtuple('TE_cluster', 'chr minpos maxpos minlength maxlength variant_ids variants TE_singles')
    inputlist=inputl.split(',')

    for inputfile in inputlist:
        i=0
        vcf = pysam.VariantFile(inputfile)

        for v in vcf:
            i=i+1
            TE_tuples = get_info(v, TE_single, TE_cluster)
            list_of_singles.append(TE_tuples[0])
            list_of_clusters.append(TE_tuples[1])
            print(v.id)

    list_of_singles=sorted(list_of_singles, key=attrgetter('pos'))
    list_of_singles=sorted(list_of_singles, key=attrgetter('chr'))

    list_of_clusters=sorted(list_of_clusters, key=attrgetter('minpos'))
    list_of_clusters=sorted(list_of_clusters, key=attrgetter('chr'))

    change=0
    old_clusters=[]
    in_clusters=[]
    new_clusters=[]
    counter = 0
    while counter < len(list_of_singles)-1:
        compare_p=list_of_clusters[counter]
        if list_of_singles[counter].variant_id in in_clusters:
            counter = counter+1
            continue
        new_cluster = compare_p
        inside_counter=1
        while True:
            if comparenear(list_of_singles[counter+inside_counter], compare_p):
                if compare(list_of_singles[counter+inside_counter], compare_p):
                    new_cluster = make_cluster(list_of_singles[counter+inside_counter], compare_p, TE_cluster)
                    compare_p=new_cluster
                inside_counter = inside_counter+1
                if inside_counter+counter >= len(list_of_singles):
                    break
            else:
                break
        counter = counter+1
        new_clusters.append(new_cluster)
        in_clusters.extend(new_cluster.variant_ids)


    output=pysam.VariantFile(outputfile, "w", header=vcf.header)
    new_info=open("merged_ids.tsv", "w+")

    soma=1
    germl=1
    for best in new_clusters:
        merged, germl, merged_info = merge_output(best, germl)
        output.write(merged)
        new_info.write(merged_info[0] + "\t" + merged_info[1] + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
