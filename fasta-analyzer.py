#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import sys
from matplotlib import pyplot as plt


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

parser.add_argument("-p", "--plot", 
                    help = "Plot GC content against length.", 
                    action = "store_true")

args = parser.parse_args()




# Each sequence becomes one object of the class Fasta.
class Fasta(object):

    def __init__(self):
        pass

    # Read fasta file.
    def read_file(self, filename):
        filename.seek(0)
        self.name, self.seq = None, []
        for line in filename:
            if line.startswith('>'):
                if self.name:
                    yield (self.name, ''.join(self.seq))
                self.name, self.seq = line, []
            else:
                self.seq.append(line)
        if self.name: yield (self.name, ''.join(self.seq))

    def header(self):
        return self.name[1:].rstrip()

    def sequence(self):
        return self.seq

    def length(self):
        self.seq = str(self.seq)
        return len(self.seq)

    # Calculates the number of G and C relative to the 
    # total number of G, C, A and T.
    def gccount(self):
        self.seq = str(self.seq)
        self.gc = self.seq.count("G") + self.seq.count("C")
        self.total = self.seq.count("G") + self.seq.count("C") + \
                     self.seq.count("A") + self.seq.count("T")
        self.content = (float(self.gc) / self.total) * 100
        return self.content

    # Print the outputs that were chosen and save the gc content 
    # and length in lists.
    def output(self):
        dictionary = {}
        for line in self.read_file(args.infile):
            dict2 = {self.header() : self}
            dictionary.update(dict2)
            if args.header == True:
                print self.header(), '\t',
            if args.length == True:
                print dictionary[self.header()].length(), '\t',
            if args.gccontent == True:
                print dictionary[self.header()].gccount(),
            if len(sys.argv) > 2:    # If no flags are given, 
                                     # no line breaks are printed.
                print	# Just there to introduce a line break.
            if args.plot == True:
                plt.scatter(dictionary[self.header()].length(), 
                dictionary[self.header()].gccount())
                plt.ylabel('GC content (%)')
                plt.xlabel('Length (nt)')
        plt.show()
        args.infile.close()




if __name__ == "__main__":
    Fasta().output()

