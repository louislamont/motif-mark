#!/usr/bin/env python3

# Motif Mark
# How is the splicing decision made?
# ISS - intron splicing silencer
# ISE - intron splicing enhancer
# ESS - exon splicing silencer
# ESE - exon splicing enhancer

# MBNL is both an activator and a repressor of exon inclusion for pre-mrna splicing
# When the MBNL binding site is downstream of the exon, it tends to cause inclusion
# When the MBNL binding site is upstream of teh exon, it tends to favor exon skipping
# RBFOX1 has a similar pattern of regulation
# MBNL and RBFOX1 have similar expression profiles
# Increased expression leads to alternative splicing changes critical for development
# Depletions in expression have been associated with neurological disorders
# Some MBNL regulated splicing events have RBFOX binding sites

# Research question: How do MBNL1 and RBFOX1 work together to regulate
# splicing events in transcripts that are regulated by both?

# Our goal: develop a python script to plot protein binding motis on an image
# of the exon and flanking sequences

# Discuss general strategy for the algorithm
# What should the input look like?
# Wht should the output look like?
# What functions are you going to write?
# How can you handle multiple motifs and multiple genes?
# How are you going to identify motifs with ambiguity, i.e. YCGY

# Output:
# Vector based image
# to scale!
# Introns vs. Exons
# Denote multiple motifs
# Need a key

# pycairo for visualization

# Functions:
# Parse FASTA
# Parse file with motifs
# Plotting function
# 

import argparse
import re
import cairo
import random

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
def draw_gene(seq, motifs):
    # Find length of seq, where introns and exons are 
    # find boundaries between introns and exons
    seqlength=len(seq)
    introns = []
    for match in re.finditer("[a-z]+", seq):
        introns.append(match.span())
        print(introns)
    exons = []
    for match in re.finditer("[A-Z]+", seq):
        exons.append(match.span())
        #print(exons)
    # Find motifs in seq, convert motifs to relative position (0-1)
    matches={}
    for onemotif in motifs:
        for match in re.finditer(motifs[onemotif], seq.upper()):
            # Add match position to dict with ambiguous motif as key
            matches.setdefault(onemotif, []).append(match.span())
            #print(match.span())
    #print(matches)
    
    WIDTH, HEIGHT = seqlength, 400

    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
    ctx = cairo.Context(surface)

    ctx.scale(WIDTH, HEIGHT)

    # Paint white background
    ctx.set_source_rgb(1, 1, 1)
    ctx.rectangle(0, 0, 1, 1)
    ctx.fill()

    # Set color back to black
    ctx.set_source_rgb(0, 0, 0)
    # Set intron line width and draw introns
    ctx.set_line_width(0.02)
    for loc in introns:
        print(loc[0], loc[1])
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
        ctx.set_font_size(0.05)
        ctx.move_to(0.1, labelloc)
        ctx.show_text(motif)
        labelloc=labelloc+0.05

        # Draw motifs
        for loc in matches[motif]:
            ctx.move_to(loc[0]/seqlength, 0.5)
            ctx.line_to(loc[1]/seqlength, 0.5)
            ctx.stroke()

    surface.write_to_png("example.png")
    surface.finish()
    # Initialize surface (300x600?) and scale
    # Draw lines for individ boundaries
    # overlay boxes for motifs on top
    # add labels and colors for motifs
    return None


iupac_dict = {"R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]", "K": "[GT]",
              "M": "[AC]", "B": "[CGT]", "D": "[AGT]", "H":"[ACT]", "V":"[ACG]",
              "N":"[ACTG]"}


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

print(iupac_motifs)

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
                print("Make first graphics thing")
                print(header, buffer)
                # Do make grapics thing (probably a function)

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
    print("We'll make final graphics here")
    print(header, buffer)
    draw_gene(buffer, iupac_motifs)
