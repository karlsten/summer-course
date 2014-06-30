#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys



# This code is partly based on https://github.com/mtop/ngs/blob/master/fp.py by Mats Töpel



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




# Each sequence becomes one object of the class Fasta.
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

    # Calculates the number of G and C relative to the 
    # total number of G, C, A and T.
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





# Print the outputs that were chosen and saves the gc content 
# and length in lists.
def output():
    for name, seq in read_file(args.infile):
        fs = Fasta(name, seq)
        if args.header == True:
            print fs.header(), '\t',
        if args.length == True:
            lengthlist.append(fs.length())
            print fs.length(), '\t',
        if args.gccontent == True:
            gclist.append(fs.gccount() * 100)
            print round(fs.gccount() * 100, 4), '%',
        if len(sys.argv) > 2:    # If no flags are given, 
                                 # no line breaks are printed.
            print	# Just there to introduce a line break.
    args.infile.close()





gclist = []
lengthlist = []





if __name__ == "__main__":
    output()

