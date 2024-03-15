#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Script merges TEs in a vcf file

import os
import sys, getopt
from typing import Sequence
import pysam
import subprocess
from functools import lru_cache
from collections import defaultdict
from cyvcf2 import VCF
from pathlib import Path
from collections import namedtuple
import mappy
import logging as log
from operator import attrgetter
from operator import itemgetter

def get_te_class(te_type):	#return TE class for different TE subtypes
    te_class=''
    if te_type.startswith("Alu"):
        te_class="Alu"
    elif te_type.startswith("L1"):
        te_class="L1"
    elif te_type.startswith("SVA"):
        te_class="SVA"
    elif te_type.startswith("HERV") or  te_type.startswith("LTR"):
        te_class="LTR"
    elif te_type.startswith("MER"):
        te_class="MER"
    elif te_type.startswith("THE"):
        te_class="THE"
    elif te_type.startswith("Tigger"):
        te_class="Tigger"
    elif te_type.startswith("MLT"):
        te_class="MLT"
    elif te_type.startswith("PRIMA4"):
        te_class="PRIMA4"
    elif te_type.startswith("MST"):
        te_class="MST"
    elif te_type.startswith("Harlequin"):
        te_class="Harlequin"
    elif te_type.startswith("LOR1"):
        te_class="LOR"
    elif te_type.startswith("PABL"):
        te_class="PABL"
    else:
        te_class=te_type
    return(te_class)

def retros(sequence, retrofasta):	#identify TE sequence
    a = mappy.Aligner(retrofasta,preset="map-ont", k=11, w=6)
    if not a:
        raise Exception("ERROR: failed to load/build index")
        sys.exit("ERROR: failed to load/build index")
    mapq=0
    besthit=None
    for hit in a.map(sequence):
        if hit.mapq >= mapq:
            besthit=hit
            mapq=hit.mapq
    te_start=besthit.r_st
    te_end=besthit.r_en
    strand=besthit.strand
    te_type=besthit.ctg
    te_class=get_te_class(te_type)
    return(te_class, te_start, te_end, strand, besthit.mapq)


def get_info(v, TE_single, TE_cluster, retrodata):	#returns info about individual insertion and insertion cluster
    list_of_TE_info = retros(v.alts[0], retrodata)
    single = TE_single(v.chrom, v.pos, list_of_TE_info[0], list_of_TE_info[2] - list_of_TE_info[1], list_of_TE_info[3], v.id, v, list_of_TE_info[4])
    cluster = TE_cluster(v.chrom, v.pos, v.pos, list_of_TE_info[0], list_of_TE_info[2] - list_of_TE_info[1], list_of_TE_info[2] - list_of_TE_info[1], list_of_TE_info[3], [v.id], [v], [single])
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
    if a.type == b.type and a.strand == b.strand and matchnumbers(a.length,b.minlength,b.maxlength, 0.75, "length") and a.variant_id not in b.variant_ids:
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
    return(TE_cluster(a.chr, minpos, maxpos, a.type, minlength, maxlength, a.strand, te_id_list, te_list, single_list))

def merge_output(cluster, soma, germl):	#create output for merged variant cluster
    id=""
    if len(cluster.variant_ids) > 1:
        id="G_"+str(germl)
        germl=germl+1
    elif not cluster.variant_ids[0].startswith("TE_"):
        id="S_"+str(soma)
        soma=soma+1
    else:
        id="G_"+str(germl)
        germl=germl+1
    sum=0
    i=0
    TES=''
    bestmapq=-1
    best_seq=''
    variant=cluster.TE_singles[0].variant.copy()
    NU=variant.info['NU']
    NG=variant.info['NG']
    NA=variant.info['NA']
    S=variant.info['S']
    TG=variant.info['TG']
    TA=variant.info['TA']

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
        if a.mapq > bestmapq:
            best_seq=v.alts[0]
            bestmapq=a.mapq
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

    variant.info.__setitem__('TE_TYPE', cluster.type)
    variant.info.__setitem__('TE_MAPQ', str(bestmapq))
    return(variant, soma, germl, [id, TES[:-1]])


def main(argv):
    outputfile=''
    inputl=''
    try:
        opts, args = getopt.getopt(argv, "hi:o:r:")
    except getopt.GetoptError:
        print("insertion_annotation.py -i <inputvcf> -o <output> -r <fasta of repeats>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("test.py -c <cram> -i <inputputlist> -r <fasta of repeats>")
            sys.exit()
        elif opt in ("-i"):
            inputl = arg
        elif opt in ("-o"):
            outputfile = arg
        elif opt in ("-r"):
            retrodata = arg


    list_of_singles=[]
    list_of_clusters=[]
    TE_single = namedtuple('TE_single', 'chr pos type length strand variant_id variant mapq')
    TE_cluster = namedtuple('TE_cluster', 'chr minpos maxpos type minlength maxlength strand variant_ids variants TE_singles')
    inputlist=inputl.split(',')

    for inputfile in inputlist:
        i=0
        vcf = pysam.VariantFile(inputfile)
        vcf.header.add_meta('INFO', items=[('ID',"TE_TYPE"), ('Number',1), ('Type','String'), ('Description','Type of transposon')])
        vcf.header.add_meta('INFO', items=[('ID',"TE_MAPQ"), ('Number',1), ('Type','String'), ('Description','Transposon mapq')])

        for v in vcf:
            i=i+1
            TE_tuples = get_info(v, TE_single, TE_cluster, retrodata)
            list_of_singles.append(TE_tuples[0])
            list_of_clusters.append(TE_tuples[1])

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
        merged, soma, germl, merged_info = merge_output(best, soma, germl)
        output.write(merged)
        new_info.write(merged_info[0] + "\t" + merged_info[1] + "\n")

if __name__ == "__main__":
    main(sys.argv[1:])

