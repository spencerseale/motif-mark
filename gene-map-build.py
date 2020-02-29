#!/usr/bin/env python

#run from conda activate snake
import argparse
import fxns
import cairo

#take in user input
def get_args():
    parser = argparse.ArgumentParser(description="Build an image of motif and genes")
    parser.add_argument("-f", "--file", help="Input fasta file", type=str, required=False)
    parser.add_argument("-m", "--motif", help="File containing motif sequences to be searched", type=str, required=False)
    return parser.parse_args()

args = get_args()

#setting input as user defined files
in_file = args.file
mot_file = args.motif

#setting surfaces
surface = cairo.PDFSurface("motif_img.pdf", 1000, 500)

#convert multi-line fasta into single line
out_fa = fxns.fasta_converter(in_file)

#get list of motifs from input motif txt file
mot_holder = fxns.motif_list(mot_file)

#create picture of exons, introns, and motifs
fxns.build_picture(surface, out_fa, mot_holder)

#building and locating legend
locator = 30
for idx, mo in enumerate(mot_holder):
    col1 = (0.7 + idx/11)
    col2 = (0.3 + idx/11)
    col3 = (0.9 + idx/11)
    fxns.draw_legend(surface, mo, 910, locator, 950, locator, col1, col2, col3)
    locator += 20

#closing the surface
surface.finish()
