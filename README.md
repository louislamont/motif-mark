# motif-mark

Motif Mark is a tool for locating and visualizing user-specified motifs within sequences. `Wright Motif Mark.py` takes a FASTA file with sequences of interest and a motif file containing one motif per line (no commas), and outputs an .svg file of the motif locations for each sequence. Example files are included. Motif Mark is capable of handling ambiguous motifs.

Usage:

```./Wright\ Motif\ Mark.py -f insr.fasta -m test.txt```

I know you mentioned you had some problems with the shebangs in the last assignment. If you encounter any trouble here, it's probably because I'm calling python3 in my shebang (since "python" for me loads python 2).