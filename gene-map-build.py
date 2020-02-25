#!/usr/bin/env python

#run from conda activate snake
import fxns
import cairo

in_file = "/Users/spencerseale/motif-mark/Figure_1.fasta"
mot_file = "/Users/spencerseale/motif-mark/Fig_1_motifs.txt"

#test motif file
#mot_file = "./mot_test.txt"

#setting surfaces
surface = cairo.PDFSurface("motif_img.pdf", 1000, 500)

#convert multi-line fasta into single line
out_fa = fxns.fasta_converter(in_file)

#get list of motifs from input motif txt file
mot_holder = fxns.motif_list(mot_file)

#edit motifs
edited_holder = fxns.edit_list(mot_holder)
#print(edited_holder)

#create picture of exons, introns, and motifs
fxns.build_picture(surface, out_fa, edited_holder)

locator = 30
for idx, mo in enumerate(edited_holder):
    col1 = (0.7 + idx/11)
    col2 = (0.3 + idx/11)
    col3 = (0.9 + idx/11)
    fxns.draw_legend(surface, mo, 910, locator, 950, locator, col1, col2, col3)
    locator += 20

surface.finish()
