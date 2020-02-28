#!/usr/bin/env python

import re
import cairo
import itertools as it

#dict holding possible ambiguities, add as many as you'd like checked
ambig = {
    "y" : ["c", "t"]
}

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

#pull motifs from input file and create list, turn rna into dna
def motif_list(motif):
    with open(motif, "r") as opmot:
        mot_holder = []
        for line in opmot:
            #i changed this
            line = line.strip()
            #replacing any uracils for thymines, only want fastas to be U
            mot_holder.append(line.replace('u', 't').replace('U', 'T'))
    return mot_holder

#getting all possible combinations for abmiguities
def combo_picker(mo, dict):
    new_list = []
    for nt in mo:
        nested_list = []
        if nt in dict:
            nested_list.append(dict[nt])
            new_list.append(nested_list)
        else:
            nested_list.append(list(nt))
            new_list.append(nested_list)
    x = list(it.product(*new_list))
    final_list = []
    for i in x[0]:
        final_list.append(i)
    poss_combo = list(it.product(*final_list))
    return ["".join(i) for i in poss_combo]

#draw motifs given their coordinates
def draw_motif(srf, motif, move_x, move_y, line_x, line_y, col1, col2, col3):
    ctx = cairo.Context(srf)
    ctx.set_line_width(len(motif))
    ctx.move_to(move_x + (len(motif)/2), move_y)
    ctx.line_to(line_x + (len(motif)/2), line_y)
    ctx.set_source_rgb(col1, col2, col3)
    ctx.stroke()

#draw legend
def draw_legend(srf, motif, leg1x, leg1y, leg2x, leg2y, col1, col2, col3):
    ctx = cairo.Context(srf)
    ctx.set_line_width(len(motif))
    ctx.set_font_size(9)
    ctx.move_to((leg1x - 80), leg1y + 3)
    ctx.show_text(motif)
    ctx.move_to(leg1x, leg1y)
    ctx.line_to(leg2x, leg2y)
    ctx.set_source_rgb(col1, col2, col3)
    ctx.stroke()

#create picture of introns, exons, and motifs in a given fasta sequence
def build_picture(srf, fasta, mot_holder):
    with open(fasta, "r") as opfasta:
        #width then height
        context = cairo.Context(srf)
        context.set_line_width(1)
        point_x = 10
        point_y = 50
        for line in opfasta:
            line = line.strip()
            if ">" not in line:
                #drawing line
                context.move_to(point_x, point_y)
                context.line_to((len(line)+point_x), point_y)
                context.set_font_size(7)
                context.move_to(point_x, point_y-13)
                context.show_text(header)
                #adding exon
                exon = re.findall(r'[A-Z]+', line)
                #find location of where exon begins, adding point_x as this is the indent
                exon_start = line.find(exon[0]) + 1 + point_x
                #drawing exon as a rectangle
                context.rectangle(exon_start, point_y-10, len(exon[0]), 20)
                #adding motifs from list
                for idx, mo in enumerate(mot_holder):
                    col1 = (0.7 + idx/11)
                    col2 = (0.3 + idx/11)
                    col3 = (0.9 + idx/11)
                    poss_mot = combo_picker(mo.lower(), ambig)
                    for x in poss_mot:
                        mot = re.findall(x, line.lower())
                    #old_slice acts to track the positioning on the line when pieces of it get removed
                    #also setting to 1 since python is 0-counting
                    #sliding by 1 position each time to find overlapping motifs
                        old_slice = 1
                        line1 = line.lower()
                        for x in range(len(mot)):
                            m_loc = line1.find(mot[x]) + old_slice
                            draw_motif(srf, mo, (m_loc+point_x), point_y-7, (m_loc+point_x), point_y+7, col1, col2, col3)
                            #string slicing away the old portion so .find() works, it can only find first occurance
                            old_slice = m_loc + 1
                            #slicing line to begin on index after where last motif ended
                            line1 = line1[(line1.find(mot[x])+1):]
                context.stroke()
                point_y += 50
            else:
                header = line
