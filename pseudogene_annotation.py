#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys, getopt
from typing import Sequence
import pysam
import subprocess
from functools import lru_cache
from collections import defaultdict
#from cyvcf2 import VCF
from pathlib import Path
from collections import namedtuple
import mappy
import logging as log
from operator import attrgetter
from operator import itemgetter
import tempfile

def read_is_target(read, coordinate):
    q_start = read.reference_start
    q_chr = read.reference_name
    c_chr = coordinate.split(':')[0]
    c_start = int(coordinate.split(':')[1])
    if q_chr == c_chr and (abs(c_start - q_start) < 10000 or abs(q_start - c_start) < 10000):
        return True
    else:
        return False

def tsd_rec(targets):
    first=targets[0]
    second=targets[0]
    tsdtype=''
    best_target=first
    sorted_targets = sorted(targets, key=lambda x: x.mapq, reverse=True)
    best_target=sorted_targets[0]
    second_best=sorted_targets[1]
    for a in [best_target, second_best]:
        if a.reference_start <= first.reference_start:
            first=a
        else:
            second=a
    if first.reference_end > second.reference_end:
        return None
    size=first.reference_end-second.reference_start
    if size > 0:
        tsdtype="target_site_duplication"
        contig1=tsd_in_contig(second, 0, size)
        contig2=tsd_in_contig(first, len(first.query_alignment_sequence)-size, len(first.query_alignment_sequence))
        list_contig=[contig1, contig2]
        coordinates = first.reference_name + ":"  + str(second.reference_start) + "-" + str(first.reference_end)
        if contig1 is None or contig2 is None:
            return None
    elif size < 0:
        tsdtype="target_site_deletion"
        contig1=tsd_in_contig(second, 0, 0)
#        contig2=tsd_in_contig(first, len(first.query_alignment_sequence)-size, ))
        list_contig=[contig1]
        coordinates = first.reference_name + ":"  + str(first.reference_end) + "-" + str(second.reference_start)
        if contig1 is None:
            return None
    if size == 0:
        return None
    return coordinates, list_contig, tsdtype

def tsd_in_contig(read, start, end):
    count = 0
    clipped = 0
    for a in read.cigartuples:
        if a[0] == 4 or a[0] == 5:
            clipped = clipped + a[1]
        if a[0] == 0 or a[0] == 1:
            count = count + a[1]
        if count >= end:
            return(clipped + start, clipped + end)

def recognize_parts(samfile, coordinates):
    target=[]
    insertion=[]
    for read in samfile:
        if read_is_target(read, coordinates):
            target.append(read)
        elif read.query is not None:
            insertion.append(read)
    return target, insertion


def recognize_polys(fasta, target_breakpoints):
    all_polys = []
    for a in fasta.references:
        sequence = fasta.fetch(a)
        letters = ("a","t")
        for letter in letters:   
            polys=is_poly(sequence, letter)
            for poly in polys:
                if near_breakpoint(poly, target_breakpoints, 10):
                    all_polys.append(poly)
        return(all_polys)

def EN_cutsite_rec(fasta, encs, annotation):
    encs_list=[]
    insertions=[]
    limit=10
    minus=1
    for a in fasta.references:
        while True:
            try:
                sequence = fasta.fetch(a, 0, fasta.lengths[0]-minus)
                break
            except ValueError:
                minus=minus+1
                
#        sequence = fasta.fetch(a, 0, fasta.lengths[0]-533)
    for a in annotation:
        if a.type == "insertion":
            insertions.append([a.start, a.end])
    for ins in insertions:
        firstseq=sequence[ins[0]-limit:ins[0]]
        secondseq=sequence[ins[1]:ins[1]+limit]
        for subsec in encs:
            if subsec[0] in firstseq and subsec[1] in secondseq:
                first_score=limit-firstseq.index(subsec[0])
                second_score=secondseq.index(subsec[1])
                encs_list.append([subsec[0], subsec[1], first_score+ins[0]-limit, first_score+len(subsec[0])+ins[0]-limit, second_score+ins[1], second_score+len(subsec[1])+ins[1], first_score+second_score])
    largest=0
    largest_en=''
    for item in encs_list:
        if item[6] > largest:
            largest=item[6]
            largest_en=item
    return largest_en

#    return(encs)

def encuts(sequence, encs_fasta):
    list_ens=[]
    for a in encs_fasta:
        list_ens.append([sequence.index(a), sequence.index(a)+len(a), a])
    return list_ens

def retros(sequence, retrofasta, annot):
    a = mappy.Aligner(sequence)
    insertions = []
    for part in annot:
        if part.type == "insertion":
            insertions.append([part.start, part.end])
    if not a:
        raise Exception("ERROR: failed to load/build index")
        sys.exit("ERROR: failed to load/build index")
    for name, seq, qual in mappy.fastx_read(retrofasta):
        for hit in a.map(seq):
        
            for ins in insertions:
                if ((hit.r_st - ins[0] > -100 and hit.r_st < ins[1]) and (hit.r_en > ins[0] and hit.r_en - ins[1] < 100)) or ((hit.r_st - ins[1] > -0 and hit.r_st < ins[0]) and (hit.r_en > ins[1] and hit.r_en -ins[0] < 100)):
                    if hit.strand == +1:
                        strand="+"
                    elif hit.strand == -1:
                        strand="-"
                    else:
                        starnd="."
                    return(name, hit.r_st, hit.r_en, hit.q_st+1, hit.q_en+1, strand)


def pseudos(sequence, pseudofasta, annot):
    insertions = []
    for part in annot:
        if part.type == "insertion":
            insertions.append([part.start, part.end])
    pseudos=[]
    Pseudo = namedtuple('Pseudo', 'transcript start_r end_r start_c end_c strand mapq alignment_score')
    output = tempfile.NamedTemporaryFile()
    os.system("minimap2 " + sequence + " " + pseudofasta + " > " + output.name)
    samput = tempfile.NamedTemporaryFile()
    os.system("minimap2 " + sequence + " " + pseudofasta + " -a > " + samput.name)
#    os.system("minimap2 " + sequence + " " + pseudofasta + " -a")

    f = open(output.name,'r')
    lines = f.readlines()
    for line in lines:
        #print(line)
        splits=line.split('\t')
        q_st=int(splits[2])
        q_en=int(splits[3])
        r_st=int(splits[7])
        r_en=int(splits[8])
        strand=splits[4]
        mapq=int(splits[11])
        alignment_score=" "
        for ins in insertions:
           if ((r_st - ins[0] > -100 and r_st < ins[1]) and (r_en > ins[0] and r_en - ins[1] < 100)): 
#or ((r_st - ins[1] > -0 and r_st < ins[0]) and (r_en > ins[1] and r_en -ins[0] < 100)):
              print(splits[0])
              asout = (subprocess.check_output("grep " + splits[0] + " " + samput.name + " | tail -n1 | cut -f14 | cut -d':' -f3", shell=True)).strip()
              if asout==b'':
                  asout=0
              pseudos.append(Pseudo(splits[0], r_st, r_en, q_st+1, q_en+1, strand, mapq, int(asout)))
#    samput = tempfile.NamedTemporaryFile()
#    os.system("minimap2 " + pseudofasta + " " + sequence + " -a > " + samput.name)
    print(pseudos)

    output.close()
    bestpseudo = select_pseudo(pseudos)
    return(bestpseudo)
#    a = mappy.Aligner(sequence)
#    insertions = []
#    for part in annot:
#        if part.type == "insertion":
#            insertions.append([part.start, part.end])
#    if not a:
#        raise Exception("ERROR: failed to load/build index")
#        sys.exit("ERROR: failed to load/build index")



#    for name, seq, qual in mappy.fastx_read(pseudofasta):
#        for hit in a.map(seq):
  #  for hit in str(minioutput):
  #      print(hit)
            #print(name)
#            print(hit)
 #           for ins in insertions:
                #print(ins)
                #print(hit.r_st)
                #print(hit.r_en)
#                if ((hit.r_st - ins[0] > -100 and hit.r_st < ins[1]) and (hit.r_en > ins[0] and hit.r_en - ins[1] < 100)) or ((hit.r_st - ins[1] > -0 and hit.r_st < ins[0]) and (hit.r_en > ins[1] and hit.r_en -ins[0] < 100)):
#                    if hit.strand == 1:
#                        strand="+"
#                    elif hit.strand == -1:
#                        strand="-"
#                    else:
#                        strand="."
#            #        print("matches")
#                    pseudos.append(Pseudo(name, hit.r_st, hit.r_en, hit.q_st+1, hit.q_en+1, strand, hit.mapq, hit.is_primary))
#    bestpseudo = select_pseudo(pseudos)
#    return(bestpseudo)

def select_pseudo(pseudos):
    bmapq=0
    #print(pseudos)
    for p in pseudos:
    #    print(p)
        if p.mapq > bmapq:
            bmapq = p.mapq
    bscore=0
    best=''
    #print(bmapq)
    for p in pseudos:
#        print("JOU")
 #       print(p)
        if p.mapq == bmapq:
#            print(p)
    #        print(p.end_r-p.start_r)
            if p.alignment_score > bscore:
                bscore = p.alignment_score
            #    blength = p.end_r-p.start_r
    #        if p.is_primary:
                best = p
        #print(blength)
    return(best)

 


def near_breakpoint(object, target_breakpoints, limit):
    near = False
    for part in target_breakpoints:
        if abs(part - object[0]) <  limit or abs(part - object[1]) < limit:
            near = True
    return near
            

def is_poly(seq, tail_letter):
    nullpos = 0
    maxval = 0
    maxpos=0
    score=0
    tails = []
    for index, ch in enumerate(seq):
        if ch.lower() == tail_letter:
            score = max([score + 1, 0])
        else:
            score = max([score - 3, 0])
        if score==0 or index == len(seq) - 1:
            if maxval > 5:
                tails.append([nullpos+1,maxpos+1,tail_letter, seq[nullpos+1:maxpos+1]])
            maxval=0
            nullpos=index
        if score >= maxval:
            maxval=score
            maxpos=index
    return tails

def target_read_insertions(target):
    count=0
    first_target=0
    second_target=0
    first_ref_end=0
    ref_count=target.reference_start
    for tuplet in target.cigartuples:
        if tuplet[0] == 1 and tuplet[1] > 20:
            first_target=count
            count = count + tuplet[1]
            second_target = count
            first_ref_end=ref_count
        elif tuplet[0] == 0:
            count = count + tuplet[1]
            ref_count = ref_count + tuplet[1]
        elif tuplet[0] == 1:
            count = count + tuplet[1]
        elif tuplet[0] == 2:
            ref_count = ref_count + tuplet[1]

    target1=[0,first_target,"target", target.reference_name + ":" + str(target.reference_start) + "-" + str(first_ref_end)]
    target2=[second_target,count,"target", target.reference_name + ":" + str(first_ref_end) + "-" + str(target.reference_end)]
    insertion=[first_target,second_target,"insertion", "."]
    return target1, insertion, target2

def find_insertions(annotation_of_contig, length_of_c):
    prevend=1

    annotation_insertions=[]
    for a in annotation_of_contig:
        if a.type == "target":
            if a.start - prevend > 5:
                annotation_insertions.append([prevend + 1, a.start - 1])
            prevend=a.end
    if length_of_c - prevend > 100:
        annotation_insertions.append([prevend + 1, length_of_c])        
    return annotation_insertions

def get_sequence(coordinates, name, fasta, read_fasta):
    sequence = (fasta.fetch(name, coordinates[0], coordinates[1]))
    read_fasta.write(">" + name + "" + str(coordinates[0]) + "" + str(coordinates[1]) + "\n")
    read_fasta.write(sequence + "\n")


def realign(read, parts, fasta, reference, breakpoint, outputdir):
    first_part = [read.query_alignment_start, parts[1]+read.query_alignment_start]
    second_part = [parts[0]+read.query_alignment_start, read.query_alignment_length+read.query_alignment_start]
    contig_parts=[first_part, second_part]
    read_fasta=open(outputdir + read.query_name +".fasta", "w+")
    for part in contig_parts:
        get_sequence(part, read.query_name, fasta, read_fasta)
    read_fasta.close()
    a = mappy.Aligner(outputdir + read.query_name +".fasta")
    if not a:
        raise Exception("ERROR: failed to load/build index")
#        sys.exit("ERROR: failed to load/build index")
    chr = breakpoint.split(':')[0]
    bp = int(breakpoint.split(':')[1])
    hits=[]
    for name, seq, qual in mappy.fastx_read(reference):
        if name == chr:
            for hit in a.map(seq): # traverse alignments
                if hit.is_primary and (abs(hit.q_st - bp) < 3000 or abs(hit.q_en - bp) < 3000 ):
                    if hit.r_st != second_part[0] and hit.r_st != first_part[0]:
                        hits.append([hit.r_st+second_part[0], hit.r_en+second_part[0], chr + ":" + str(hit.q_st) + "-" + str(hit.q_en), hit.mapq, hit])
                    else:
                        hits.append([hit.r_st, hit.r_en, chr + ":" + str(hit.q_st) + "-" + str(hit.q_en), hit.mapq, hit])
    hit0=hits[0]
    first_hit=hit0[4]
    second_hit=hit0[4]
    for hit1 in hits:
        if hit1[4].q_st <= first_hit.q_st:
            first_hit=hit1[4]
        else:
            second_hit=hit1[4]
    size=first_hit.q_en-second_hit.q_st
    if size > 0:        
        tsdtype="target_site_duplication"
        contig1=tsd_in_contig_mappy(second_hit.cigar, second_hit.r_st+second_part[0], second_hit.r_st+size+second_part[0], second_part[0])
        contig2=tsd_in_contig_mappy(first_hit.cigar, first_hit.r_en-size, first_hit.r_en, second_part[0])
        list_contig=[contig1, contig2]
        coordinates = chr + ":"  + str(second_hit.q_st) + "-" + str(first_hit.q_en)
        if contig1 is None or contig2 is None:
            return None
    elif size < 0:
        tsdtype="target_site_deletion"
        #contig1=tsd_in_contig_mappy(second_hit.cigar, 0, 0, 0)
        contig1=(parts[0]+read.query_alignment_start+1, parts[0]+read.query_alignment_start+1)
        list_contig=[contig1]
        coordinates = chr + ":"  + str(first_hit.q_en) + "-" + str(second_hit.q_st)
    else:
        return None

    return hits, [coordinates, list_contig, tsdtype]

def tsd_in_contig_mappy(cigar, start, end, previous):

    count = previous
    clipped = 0
    for a in cigar:
        print(a)
        if a[1] == 4 or a[1] == 5:
            clipped = clipped + a[0]
        if a[1] == 0 or a[1] == 1:
            count = count + a[0]
        if count >= end:
            return(clipped + start, clipped + end)
    

def aligned_part(target):
    count=0
    start=0
    end=0
    for cigar in target.cigartuples:
        if cigar[0] == 4 or cigar[0] == 5:
            if count==0:
                start=cigar[1]
                count=cigar[1]
            else:
                end=count
        elif cigar[0] == 0 or cigar[0] == 1:
            count = count + cigar[1]
    if end == 0:
        end = count
    return start, end 


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hf:s:o:p:r:g:")
    except getopt.GetoptError:
        print("insertion_annotation.py [-S] -s <sam> -f <fasta> -p <position> -r <retrotransposon fasta>  -o <output>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("test.py -c <cram> -i <inputputlist>")
            sys.exit()
        elif opt in ("-s"):
            sam = arg
        elif opt in ("-f"):
            fasta = arg
        elif opt in ("-p"):
            position = arg
        elif opt in ("-r"):
            retrofasta = arg
        elif opt in ("-o"):
            outputfile = arg
        elif opt in ("-g"):
            pseudofasta = arg

    output = open(outputfile, "w+")

    reference_genome="/mnt/cg8/reference-genomes/GRCh38_no_alt/GRCh38_no_alt.fasta"
    outputdir=fasta.split('q/', 1)[0]+"q"
    fastafile = pysam.FastaFile(fasta)
    breakpoint = position
    samfile=pysam.AlignmentFile(sam, "rb")
    annotation_of_contig = []

    parts = recognize_parts(samfile, breakpoint)
    targets=parts[0]
    insertions=parts[1]
    if parts == ([],[]):

        exit()
    endonuclease_cuts=[["TTTT","A"], ["T","AAAA"], ["TTTT","G"], ["C","AAAA"], ["TTTC","A"], ["T","GAAA"]]

#    print(parts)
    
    alignments=parts[0]
    if alignments:
        read=alignments[0]
    else:
        read=parts[1][0]


    Part = namedtuple('Part', 'type contig start end content score strand')
    realigned_parts = None 


    for target in targets:
        if target.is_reverse:
           reverse="-"
        else:
            reverse="+"
        aligned_target = aligned_part(target)
        insertion_in_target = target_read_insertions(target)
        for parts in insertion_in_target:
            if parts[2] == "insertion" and parts[1] != 0:
                realigned_parts=realign(target, parts, fastafile, reference_genome, breakpoint, outputdir)            
        if insertion_in_target[0][1] != 0:
            if realigned_parts:
                for p in realigned_parts[0]:
                    annotation_of_contig.append(Part("target", target.query_name, p[0], p[1], p[2], p[3], reverse))
            else:
                for part in insertion_in_target:
                    if part[2] == "insertion":
                        annotation_of_contig.append(Part(part[2], target.query_name, part[0], part[1], part[3], ".", "."))
                    else:
                        annotation_of_contig.append(Part(part[2], target.query_name, part[0], part[1], part[3], target.mapq, reverse))            

        else:
            annotation_of_contig.append(Part("target", target.query_name, aligned_target[0], aligned_target[1], target.reference_name + ":" + str(target.reference_start) + "-" + str(target.reference_end), target.mapq, reverse))
    
    target_limits=[]
    for part in annotation_of_contig:
        if part.type == "target":
            if part.start != 0:
                target_limits.append(part.start)
            if part.end != fastafile.lengths[0]:
                target_limits.append(part.end)



    if targets:
        if insertion_in_target[0][1] != 0:
            for part2 in insertion_in_target:
                if part2[2] != "insertion":
                    if part2[0] != 0:
                       target_limits.append(part2[0])
                    if part2[1] != fastafile.lengths[0] and part2[1] != 0:
                       target_limits.append(part2[1])

    if insertions:
        for insertion in insertions:
            if insertion.is_reverse:
               reverse="-"
            else:
                reverse="+"
            aligned_ins = aligned_part(insertion)
            annotation_of_contig.append(Part("non_target", insertion.query_name, aligned_ins[0], aligned_ins[1], insertion.reference_name + ":" + str(insertion.reference_start) + "-" + str(insertion.reference_end), insertion.mapq, reverse))


    annotation_of_contig = sorted(annotation_of_contig, key=lambda x: x.end)
    annotation_of_contig = sorted(annotation_of_contig, key=lambda x: x.start)    

#    print(fastafile.lengths)
    new_insertions = find_insertions(annotation_of_contig, fastafile.lengths[0])
    for new_insertion in new_insertions:
        annotation_of_contig.append(Part("insertion", read.query_name, new_insertion[0], new_insertion[1], ".", ".", "."))
      
 
    if len(targets) > 1:
        tsd = tsd_rec(targets)
        if tsd is not None:
            for tsdpart in tsd[1]:
                annotation_of_contig.append(Part(tsd[2], read.query_name, tsdpart[0], tsdpart[1], tsd[0], ".", "."))

    if realigned_parts:
        tsd = realigned_parts[1]
        if tsd is not None:
            for tsdpart in tsd[1]:
                annotation_of_contig.append(Part(tsd[2], read.query_name, tsdpart[0], tsdpart[1], tsd[0], ".", "."))



    retrotransposons=retros(fasta, retrofasta, annotation_of_contig)

    pseudogenes=pseudos(fasta, pseudofasta, annotation_of_contig)

    if retrotransposons:
        annotation_of_contig.append(Part("transposon", read.query_name, retrotransposons[1], retrotransposons[2], retrotransposons[0] + ":" + str(retrotransposons[3]) + "-" + str(retrotransposons[4]), ".", retrotransposons[5]))

    if pseudogenes:
        print("ON")
        annotation_of_contig.append(Part("pseudogene", read.query_name, pseudogenes[1], pseudogenes[2], pseudogenes[0] + ":" + str(pseudogenes[3]) + "-" + str(pseudogenes[4]), str(pseudogenes[6]), pseudogenes[5]))


    cuts = EN_cutsite_rec(fastafile, endonuclease_cuts, annotation_of_contig)
    if cuts:
        annotation_of_contig.append(Part("endonuclease_cut_site", read.query_name, cuts[2], cuts[3], cuts[0], ".", "."))
        annotation_of_contig.append(Part("endonuclease_cut_site", read.query_name, cuts[4], cuts[5], cuts[1], ".", "."))

    polys = recognize_polys(fastafile, target_limits)
    if polys:
        for poly in polys:
            if poly[2] == "a":
                annotation_of_contig.append(Part("polyA", read.query_name, poly[0], poly[1], poly[3], ".", "."))
            if poly[2] == "t":
                annotation_of_contig.append(Part("polyT", read.query_name, poly[0], poly[1], poly[3], ".", "."))

 #       annotation_of_contig.append(Part(, read.query_name, cuts[2], cuts[3], cuts[0], ".", "."))

    annotation_of_contig = sorted(annotation_of_contig, key=lambda x: x.end)
    annotation_of_contig = sorted(annotation_of_contig, key=lambda x: x.start)

    print(annotation_of_contig)

    output.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tcontent\n")

    for annot_part in annotation_of_contig:
        output.write(annot_part.contig + "\t" + str(annot_part.start + 1) + "\t" + str(annot_part.end + 1) + "\t" + annot_part.type +  "\t" + str(annot_part.score) + "\t" + annot_part.strand + "\t" + annot_part.content + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])

