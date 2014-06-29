#!/usr/bin/env python

import argparse
import sys



# This code is partly based on https://github.com/mtop/ngs/blob/master/fp.py 



parser = argparse.ArgumentParser(description = 
                                 """Calculate length and/or GC content for 
                                    sequences from a fasta file. The result 
                                    is printed to stdout and saved in lists 
                                    'gclist' for GC content and 'lengthlist'
                                    for length""")

parser.add_argument("infile", 
                    type = argparse.FileType('r'), 
                    help = "Infile in fasta format.", 
                    default = sys.stdin)

parser.add_argument("-hd", "--header", 
                    help = "Print sequence headers.", 
                    action = "store_true")

parser.add_argument("-l", "--length", 
                    help = "Calculate the length of each sequence.",
                    action = "store_true")

parser.add_argument("-gc", "--gccontent", 
                    help = "Calculate the GC content of each sequence.", 
                    action = "store_true")

args = parser.parse_args()





class Fasta(object):

    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

    def header(self):
        return self.name[1:].rstrip()

    def sequence(self):
        return self.seq

    def length(self):
        return len(self.seq)

    # Calculates the number of G and C relative to the total number of G, C, A and T.
    def gccount(self):
        self.gc = self.seq.count("G") + self.seq.count("C")
        self.total = self.seq.count("G") + self.seq.count("C") + \
                     self.seq.count("A") + self.seq.count("T")
        self.content = float(self.gc) / self.total
        return self.content





# Read fasta file.
def read_file(infile):
    args.infile.seek(0)
    name, seq = None, []
    for line in args.infile:
        if line.startswith('>'):
            if name:
                yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))





# Adds the length of each sequence to lengthlist.
def length():
    for name, seq in read_file(args.infile):
        fs = Fasta(name, seq)
        lengthlist.append(fs.length())





# Adds the GC content of each sequence (in percent) to gclist.
def gccontent():
    for name, seq in read_file(args.infile):
        fs = Fasta(name, seq)
        gclist.append(fs.gccount() * 100)





# Print the outputs that were chosen.
def print_output():
    for name, seq in read_file(args.infile):
        fs = Fasta(name, seq)
        if args.header == True:
            print fs.header(), '\t',
        if args.length == True:
            print fs.length(), '\t',
        if args.gccontent == True:
            print round(fs.gccount() * 100, 4), '%',
        if len(sys.argv) > 2:    # If no flags are given, no line breaks are printed.
            print	# Just there to introduce a line break.





gclist = []
lengthlist = []

if args.length == True:
    length()

if args.gccontent == True:
    gccontent()

print_output()

args.infile.close()
