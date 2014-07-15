#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
from matplotlib import pyplot as plt



# This code is partly based on https://github.com/mtop/ngs/blob/master/fp.py
# by Mats Töpel



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

parser.add_argument("coverage",
                    type = argparse.FileType('r'),
                    help = "Coverage file (name and coverage separated by tab)",
                    nargs = '?')

parser.add_argument("-hd", "--header", 
                    help = "Print sequence headers.", 
                    action = "store_true")

parser.add_argument("-l", "--length", 
                    help = "Calculate the length of each sequence.",
                    action = "store_true")

parser.add_argument("-gc", "--gccontent", 
                    help = "Calculate the GC content of each sequence.", 
                    action = "store_true")

parser.add_argument("-p", "--lenplot", 
                    help = "Plot GC content against length.", 
                    action = "store_true")

parser.add_argument("-c", "--covplot",
                    help = "Plot GC content against coverage.",
                    action = "store_true")

args = parser.parse_args()





# Each sequence becomes one object of the class Fasta.
class Fasta(object):

    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        self.cov = float('nan')	# For missing coverage values, 
                                # covplot will not plot.

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
        self.content = (float(self.gc) / self.total) * 100
        return self.content

    def coverage(self):
        return self.cov





# Read coverage file and add coverage for each contig to the dictionary.
def read_covfile(infile):
    infile.seek(0)
    for line in infile:
        name = line.split('\t')[0]
        if name in dictionary.keys():
            if line.split('\t')[1].rstrip() is not '': # Avoid error caused by
                                                       # '' not being an int
                                                       # if coverage, but not 
                                                       # \t, is missing in file
                dictionary[name].cov = int(line.split('\t')[1].rstrip())
    return dictionary





# Read fasta file.
def read_file(infile):
    infile.seek(0)
    name, seq = None, []
    for line in infile:
        if line.startswith('>'):
            if name:
                fs = Fasta(name, ''.join(seq))
                dict2 = {fs.header() : fs}
                dictionary.update(dict2)
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        fs = Fasta(name, ''.join(seq))
        dict2 = {fs.header() : fs}
        dictionary.update(dict2)
    return dictionary





# Plot GC content against length. 
def lenplot():
    for key in dictionary:
        plt.scatter(dictionary[key].length(), 
        dictionary[key].gccount())
    plt.suptitle('GC - Length', fontsize = 20)
    plt.ylabel('GC content (%)')
    plt.xlabel('Length (nt)')
    plt.show()





# Plot GC content against coverage.
def covplot():
    for key in dictionary:
        plt.scatter(dictionary[key].gccount(), 
        dictionary[key].coverage())
    plt.suptitle('GC - Coverage', fontsize = 20)
    plt.ylabel('Coverage')
    plt.xlabel('GC content (%)')
    plt.show()





# Print the outputs that were chosen.
def main():
    read_file(args.infile)
    if args.coverage:
        read_covfile(args.coverage)
    for key in dictionary:
        if args.header == True:
            print dictionary[key].header(), '\t',
        if args.length == True:
            print dictionary[key].length(), '\t',
        if args.gccontent == True:
            print dictionary[key].gccount(), '\t',
        if args.coverage:
            print dictionary[key].coverage(),
        if len(sys.argv) > 2:    # If no flags are given, 
                                 # no line breaks are printed.
            print	# Just there to introduce a line break.
    if args.lenplot == True:
        lenplot()
    if args.covplot == True:
        try:
            covplot()
        except:
            print "ERROR: Correct coverage file not supplied."
    args.infile.close()





dictionary = {}





if __name__ == "__main__":
    main()

