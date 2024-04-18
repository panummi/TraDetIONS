#!/usr/bin/env python
#-*- coding: utf-8 -*-

#Script recognizes if sequence contains TE and filters out simple repeats

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Annotate called vcf inserts with repeat matching info")

    parser.add_argument("-i", "--input",
                        help="Input file [default:%(default)s]",
                        default="/dev/stdin")

    parser.add_argument("-a", "--annot",
                        help="Annotation fasta file [default:%(default)s]",
                        default=None)

    parser.add_argument("-o", "--output",
                        help="Output prefix [default:%(default)s]",
                        default="/dev/stdout")

    parser.add_argument("-V", "--verbose",default=False,const=True,nargs="?",
                        help="Be more verbose with output [and log to a file] [default:%(default)s]" )

    args = parser.parse_args()

    import logging
    if args.verbose:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s:%(funcName)s:%(levelname)s:%(message)s')
        if args.verbose!=True:
            log_file_handler = logging.FileHandler(args.verbose)
            log_file_handler.setFormatter(logging.getLogger().handlers[0].formatter)
            logging.getLogger().addHandler(log_file_handler)

    return args

class SnifflesSeqAnnotator(object):
    def __init__(self, input_file):
        from cyvcf2 import VCF
        self._vcfstrm = VCF(input_file)

    def _iter_insertions(self):
        import logging as log
        import re
        nucl_pat = re.compile("^[ACGT]+$")
        for variant in self._vcfstrm:
            seq = variant.ALT[0]
            if nucl_pat.match(seq) is None:
                log.info("Skipping {}".format(variant))
                continue

            yield variant

    def _iter_indels(self):
        import logging as log
        import re
        nucl_pat = re.compile("^[ACGT]+$")
        for variant in self._vcfstrm:
            alt = variant.ALT[0]
            ref = variant.REF
            if ref=="N" and nucl_pat.match(alt) is not None:
                s = alt
            #elif alt=="N" and nucl_pat.match(ref) is not None:
            #    s = ref
            else:
                s = ""
            yield variant,s

    def is_poly(self,seq,max_len=5):	#recognizes if sequence is a simple repeat
        from collections import Counter
        if len(seq) > 0:
            for l in range(1,max_len+1):
                base_count = Counter(seq[i:i+l] for i in range(0,len(seq),l))
                (b,c), = base_count.most_common(1)
                if len(seq)*0.9/l < c:
                        return "poly{}".format(b)
            else:
                return False

    def annotate_retros(self,retrofasta,preset="map-ont"):	#annotates transposons
        import mappy
        self.aln = mappy.Aligner(retrofasta,preset="map-ont")
        for variant in self._iter_insertions():
            unmapped=True
            RE=variant.INFO.get("RE")
            is_poly = self.is_poly(variant.ALT[0])
            if not is_poly:
                for aln in self.aln.map(variant.ALT[0]):
                    unmapped=False
                    print("{RE}\t{v.CHROM}\t{v.POS}\t{v.ID}\t{0}\t{v.ALT[0]}".format(str(aln),RE=RE,v=variant))
                if unmapped or True:
                    print(RE,str(variant).strip())
            else:
                print("{RE}\t{v.CHROM}\t{v.POS}\t{v.ID}\t{0}\t{v.ALT[0]}".format(is_poly,RE=RE,v=variant))


    def to_annotated_vcf(self,retrofasta,out_vcf,preset="map-ont",filter_fn=lambda v:True):	#produces output vcf
        import mappy
        import logging as log
        self.aln = mappy.Aligner(retrofasta,preset="map-ont", k=11, w=6)
        from cyvcf2 import VCF, Writer
        self._vcfstrm.add_info_to_header({"ID":"SVclass",
            "Description":"Class of SV (polyX, Transposable Element)",
            "Type":"String","Number":1})
        self._vcfstrm.add_info_to_header({"ID":"POSITION",
            "Description":"Position of TE in the query (middle, beginning, end, whole)",
            "Type":"String","Number":1})
        self._vcfstrm.add_info_to_header({"ID":"FLANKINGstart",
            "Description":"Number of basepairs before transposon sequence",
            "Type":"Integer","Number":1})
        self._vcfstrm.add_info_to_header({"ID":"FLANKINGend",
            "Description":"Number of basepairs after transposon sequence",
            "Type":"Integer","Number":1})
        self._vcfstrm.add_info_to_header({"ID":"poly",
            "Description":"Is element poly",
            "Type":"String","Number":1})
        self._vcfstrm.add_info_to_header({"ID":"TE_length",
            "Description":"Length of TE sequence",
            "Type":"Integer","Number":1})
        ovcf = Writer(out_vcf,self._vcfstrm)

        for variant,allele in self._iter_indels():
            unmapped=True
            RE=variant.INFO.get("RE")
            if RE is not None and RE < 3:
                continue
            is_poly = self.is_poly(allele)
            if is_poly is False:
                bestscore=int(0)
                primary=0
                for aln in self.aln.map(allele):
                   if aln.is_primary:
                      primary=1
                      if aln.mapq >= bestscore:
                        bestscore=aln.mapq
                        unmapped=False
                        SVclass=aln.ctg.split("#")[-1]
                        FLANKINGstart=aln.q_st
                        FLANKINGend=len(allele) - int(aln.q_en)
                        TE_length=aln.q_en-aln.q_st
                        if FLANKINGstart < 200 and FLANKINGend < 200:
                          POSITION="whole"
                        elif FLANKINGstart < 200:
                          POSITION="beginning"
                        elif FLANKINGend < 200:
                          POSITION="end"
                        else:
                          POSITION="middle"
                if primary == 1:
                    variant.INFO["POSITION"]=POSITION
                    variant.INFO["FLANKINGstart"]=FLANKINGstart
                    variant.INFO["FLANKINGend"]=FLANKINGend
                    variant.INFO["SVclass"]=SVclass
                    variant.INFO["TE_length"]=TE_length
                    ovcf.write_record(variant)
                variant.INFO["poly"] = str(is_poly)


if __name__ == '__main__':
    args=main()
    annot = SnifflesSeqAnnotator(args.input)
    import sys
    vcf = args.output + ".vcf"
    annot.to_annotated_vcf(args.annot, vcf)
