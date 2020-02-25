#!/usr/bin/env python

import re
import cairo

#line converter for fasta
def fasta_converter(IN_FILE):
    OUT_FILE = IN_FILE+"_ADJ.fasta"
    with open(OUT_FILE, "w") as outfile:
        with open(IN_FILE, "r") as infile:
            count = 0
            for line in infile:
                if count < 1:
                    outfile.write(line.replace('u', 't').replace('U', 'T'))
                    count += 1
                else:
                    if ">" in line:
                        outfile.write("\n")
                        outfile.write(line.replace('u', 't').replace('U', 'T'))
                    else:
                        outfile.write(line.strip().replace('u', 't').replace('U', 'T'))
    return OUT_FILE

#mine motifs from input file
def motif_list(motif):
    with open(motif, "r") as opmot:
        mot_holder = []
        for line in opmot:
            line = line.strip()
            #replacing any uracils for thymines, only want fastas to be U
            mot_holder.append(line.replace('u', 't').replace('U', 'T'))
    return mot_holder


#draw motifs given their coordinates
def draw_motif(srf, motif, move_x, move_y, line_x, line_y):
    ctx = cairo.Context(srf)
    ctx.set_line_width(len(motif))
    # ctx.move_to(move_x, move_y)
    # ctx.line_to(line_x, line_y)
    ctx.move_to(move_x+(len(motif)/2), move_y)
    ctx.line_to(line_x+(len(motif)/2), line_y)
    ctx.stroke()

#create picture of introns, exons, and motifs in a given fasta sequence
def build_picture(fasta, output, mot_holder):
    with open(fasta, "r") as opfasta:
        print(f"motifs: {mot_holder}")
        #width then height
        surface = cairo.PDFSurface(output, 1000, 500)
        context = cairo.Context(surface)
        context.set_line_width(1)
        point_x = 10
        point_y = 50
        for line in opfasta:
            line = line.strip()
            if ">" not in line:
                #drawing line
                context.move_to(point_x, point_y)
                context.line_to((len(line)+point_x), point_y)
                #adding exon
                exon = re.findall(r'[A-Z]+', line)
                #find location of where exon begins, adding point_x as this is the indent
                exon_start = line.find(exon[0]) + 1 + point_x
                #drawing exon as a rectangle
                context.rectangle(exon_start, point_y-5, len(exon[0]), 10)
                #adding motifs from list
                for mo in mot_holder:
                    mot = re.findall(mo, line)
                    print(mot)
                    #old_slice acts to track the positioning on the line when pieces of it get removed
                    #also setting to 1 since python is 0-counting
                    #sliding by 1 position each time to find overlapping motifs
                    old_slice = 1
                    line1 = line
                    for x in range(len(mot)):
                        m_loc = line1.find(mot[x]) + old_slice
                        draw_motif(surface, mo, (m_loc+point_x), point_y-9, (m_loc+point_x), point_y+9)
                        #string slicing away the old portion so .find() works, it can only find first occurance
                        old_slice = m_loc + 1
                        #slicing line to begin on index after where last motif ended
                        line1 = line1[(line1.find(mot[x])+1):]

                context.stroke()
                point_y += 50
        surface.finish()
