#!/usr/bin/env python3

import argparse
import re
import cairo
import random

# ArgParse
def get_arguments():
    #Give description of motif mark
    parser = argparse.ArgumentParser(description="Finds location of motifs in \
supplied sequences.")
    parser.add_argument("-f", "--file", help="Takes a FASTA file containing one \
or more sequences. Introns should be lower case, exons upper case.",
                        required=True, type=argparse.FileType("r"))
    parser.add_argument("-m", "--motifs", help="Takes a file containing motifs \
to search (one per line).", required=True, type=argparse.FileType("r"))
    return parser.parse_args()

args=get_arguments()

# Function to convert IUPAC characters to regex strngs
def motif_change(motifs):
    iupac_motifs={}
    for one in motifs:
        for key in iupac_dict:
            if key in one:
                replace = one.replace(key, iupac_dict[key])
        iupac_motifs[one]=replace
    return iupac_motifs

# Drawing function
def draw_gene(header, seq, motifs, counter):
    # Find length of seq and where introns and exons are
    seqlength=len(seq)
    introns = []
    for match in re.finditer("[a-z]+", seq):
        introns.append(match.span())
    exons = []
    for match in re.finditer("[A-Z]+", seq):
        exons.append(match.span())

    # Find motifs in seq, convert motifs to relative position (0-1)
    matches={}
    for onemotif in motifs:
        for match in re.finditer(motifs[onemotif], seq.upper()):
            # Add match position to dict with ambiguous motif as key
            matches.setdefault(onemotif, []).append(match.span())
    
    WIDTH, HEIGHT = seqlength, 400

    surface = cairo.SVGSurface(str(counter)+".svg", WIDTH, HEIGHT)
    ctx = cairo.Context(surface)

    ctx.scale(WIDTH, HEIGHT)

    # Paint white background
    ctx.set_source_rgb(1, 1, 1)
    ctx.rectangle(0, 0, 1, 1)
    ctx.fill()

    # Set color back to black
    ctx.set_source_rgb(0, 0, 0)

    # Add gene name
    ctx.select_font_face("Serif", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
    ctx.set_font_size(0.02)
    #ctx.FontOptions.set_width(0.01)
    ctx.move_to(0.025, 0.05)
    ctx.show_text(header)
    
    # Set intron line width and draw introns
    ctx.set_line_width(0.02)
    for loc in introns:
        ctx.move_to(loc[0]/seqlength, 0.5)
        ctx.line_to(loc[1]/seqlength, 0.5)
        ctx.stroke()

    # Same for exons
    ctx.set_line_width(0.075)
    for loc in exons:
        ctx.move_to(loc[0]/seqlength, 0.5)
        ctx.line_to(loc[1]/seqlength, 0.5)
        ctx.stroke()

    # Make location variable for motif label
    labelloc = 0.1
    # Iterate through motif dictionary
    for motif in matches:
        # Set color for motif and add label. Add small value to location
        ctx.set_source_rgb(random.uniform(0, 1), random.uniform(0, 1),
                        random.uniform(0, 1))
        ctx.select_font_face("Serif")
        ctx.set_font_size(0.025)
        ctx.move_to(0.025, labelloc)
        ctx.show_text(motif)
        labelloc=labelloc+0.035

        # Draw motifs
        for loc in matches[motif]:
            ctx.move_to(loc[0]/seqlength, 0.5)
            ctx.line_to(loc[1]/seqlength, 0.5)
            ctx.stroke()

    #surface.write_to_svg("example.svg")
    surface.finish()
    # Initialize surface (300x600?) and scale
    # Draw lines for individ boundaries
    # overlay boxes for motifs on top
    # add labels and colors for motifs
    counter=counter+1
    return counter


# Create ambiguous IUPAC dictionary
iupac_dict = {"R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]", "K": "[GT]",
              "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H":"[ACT]", "V":"[ACG]",
              "N":"[ACTG]"}

# Set counter for filenames
counter = 1

### Get Motifs
motifs=[]

# For each line in the motif file, convert to upper case
with args.motifs as motiffile:
    for line in motiffile:
        line = line.strip()
        line = line.upper()
        motifs.append(line)
# then generate regex strings for each motif
iupac_motifs=motif_change(motifs)

### Read through sequence file
with args.file as seqfile:

    firsttime=True
    header=""
    buffer=""
    # read through file line by line
    for line in seqfile:
        # If the first char of line is ">" (we're on a name line)
        if line[0] == '>':
            # If we're not on the first entry, we've got sequence in
            # the line before, so print a newline, then the next name
            if firsttime==False:
                # we're on a second+ header -> make graphic for previous buffer
                counter=draw_gene(header, buffer, iupac_motifs, counter)
                
                # Reset header and buffer
                header = line.strip()
                buffer = ""
            # Otherwise, we are on the first entry, so record header and set flag
            else:
                header = line.strip()
                firsttime=False
        # If the first char isn't a ">", we're on a sequence line
        else:
            # so add that to our buffer
            buffer = buffer + line.strip()
    # We have reached EOF. Make graphic for final buffer (function).
    counter=draw_gene(header, buffer, iupac_motifs, counter)
