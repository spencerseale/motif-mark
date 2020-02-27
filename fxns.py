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

#pull motifs from input file
def motif_list(motif):
    with open(motif, "r") as opmot:
        mot_holder = []
        for line in opmot:
            #i changed this
            line = line.strip()
            #replacing any uracils for thymines, only want fastas to be U
            mot_holder.append(line.replace('u', 't').replace('U', 'T'))
    return mot_holder

#edit motifs in list to replace those with ambiguous characters
def edit_list(motif_list):
    edited_list = []
    for element in motif_list:
        if "y" in element:
            new1 = element.replace("y", "c")
            new2 = element.replace("y", "t")
            edited_list.append(new1)
            edited_list.append(new2)
        elif "Y" in element:
            new1 = element.replace("Y", "C")
            new2 = element.replace("Y", "T")
            edited_list.append(new1)
            edited_list.append(new2)
        else:
            edited_list.append(element)
    return edited_list

#draw motifs given their coordinates
def draw_motif(srf, motif, move_x, move_y, line_x, line_y, col1, col2, col3):
    ctx = cairo.Context(srf)
    ctx.set_line_width(len(motif))
    # ctx.move_to(move_x, move_y)
    # ctx.line_to(line_x, line_y)
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
        #print(f"motifs: {mot_holder}")
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
                context.rectangle(exon_start, point_y-5, len(exon[0]), 10)
                #adding motifs from list
                for idx, mo in enumerate(mot_holder):
                    col1 = (0.7 + idx/11)
                    col2 = (0.3 + idx/11)
                    col3 = (0.9 + idx/11)
                    mot = re.findall(mo, line)

                    #line.find(mo)
                    #old_slice acts to track the positioning on the line when pieces of it get removed
                    #also setting to 1 since python is 0-counting
                    #sliding by 1 position each time to find overlapping motifs
                    old_slice = 1
                    line1 = line.lower()
                    #while mo.lower() in line1:
                    for x in range(len(mot)):
                        m_loc = line1.find(mot[x].lower()) + old_slice
                        draw_motif(srf, mo, (m_loc+point_x), point_y-9, (m_loc+point_x), point_y+9, col1, col2, col3)
                        #string slicing away the old portion so .find() works, it can only find first occurance
                        old_slice = m_loc + 1
                        #slicing line to begin on index after where last motif ended
                        line1 = line1[(line1.find(mot[x].lower())+1):]
                context.stroke()
                point_y += 50
            else:
                header = line
