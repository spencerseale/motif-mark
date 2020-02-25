#!/usr/bin/env python

#run from conda activate snake
import fxns

in_file = "/Users/spencerseale/motif-mark/Figure_1.fasta"
#mot_file = "/Users/spencerseale/motif-mark/Fig_1_motifs.txt"

#test motif file
mot_file = "./mot_test.txt"

#convert multi-line fasta into single line
out_fa = fxns.fasta_converter(in_file)

#get list of motifs from input motif txt file
mot_holder = fxns.motif_list(mot_file)

#create picture of exons, introns, and motifs
fxns.build_picture(out_fa, "test.pdf", mot_holder)
