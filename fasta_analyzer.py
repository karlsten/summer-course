#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import numpy as np
from matplotlib import pyplot as plt


# This section is based on code from 
# https://github.com/mtop/speciesgeocoder/blob/master/geocoder.py 
# by Mats Töpel:
# (also, it's not done yet)
#try:
#    import argparse
#except ImportError:
#    print >> sys.stderr, 'ERROR: The module "argparse" is not installed\n'
#    answer = raw_input('Would you like to install it now using "sudo pip \
#                        install argparse"? [y/n]: ')
#    if answer == 'y' or 'Y':
#        print 'Running "sudo pip install argparse"...'
#        from subprocess import call
#        call(["sudo", "pip", "install", "argparse"])
#    else:
#        sys.exit("ERROR: Exiting ____________________")

#try:
#    from matplotlib import pyplot as plt
#except ImportError:
#    print >> sys.stderr, 'ERROR: The library "matplotlib" is not installed\n'
#    answer = raw_input('Would you like to install it now using "sudo pip \
#                        install matplotlib"? [y/n]: ')
#    if answer == 'y' or 'Y':
#        print 'Running "sudo pip install matplotlib"...'
#        from subprocess import call
#        call(["sudo", "pip", "install", "matplotlib"])
#    else:
#        sys.exit("ERROR: Exiting ____________________")





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

parser.add_argument("-n", "--ncontent",
                    help = "Calculate the number of N in each sequence.",
                    action = "store_true")

parser.add_argument("-lg", "--lengcplot", 
                    help = "Plot GC content against length.", 
                    action = "store_true")

parser.add_argument("-cg", "--covgcplot",
                    help = "Plot GC content against coverage.",
                    action = "store_true")

parser.add_argument("-cl", "--covlenplot",
                    help = "Plot coverage against length.",
                    action = "store_true")

parser.add_argument("-ch", "--covhistogram",
                    help = " Histogram over coverage.",
                    action = "store_true")

parser.add_argument("-lh", "--lenhistogram",
                    help = "Histogram over length.",
                    action = "store_true")

args = parser.parse_args()





# Each sequence becomes one object of the class Fasta.
class Fasta(object):

    def __init__(self, name, seq):
        self.name = name[1:].rstrip()
        self.seq = seq.lower()
        self.cov = float('nan')	# For missing coverage values, 
                                # covgcplot will not plot.

    def header(self):
        return self.name

    def sequence(self):
        return self.seq

    def length(self):
        return len(self.seq)

    # Calculates the number of G and C relative to the 
    # total number of G, C, A and T.
    def gccount(self):
        self.gc = self.seq.count("g") + self.seq.count("c")
        self.total = self.seq.count("g") + self.seq.count("c") + \
                     self.seq.count("a") + self.seq.count("t")
        if self.total == 0:
            self.content = float('nan')
        else:
            self.content = round((float(self.gc) / self.total) * 100, 3)
        return self.content

    def coverage(self):
        return self.cov

    def ncontent(self):
        return self.seq.count("n")





class conname(object):
    def __init__(self, name):
        self.name = name





# Read coverage file and add coverage for each contig to the dictionary.
def read_covfile(infile, dictionary):
    infile.seek(0)
    for line in infile:
        name = line.split('\t')[0]
        if name in dictionary.keys():
            try:
                dictionary[name].cov = float(line.split('\t')[1].rstrip())
            except:
                dictionary[name].cov = float('nan')
    return dictionary





# Read fasta file.
def read_file(infile):
    dictionary = {} 
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
def lengcplot(dictionary):
    contig = conname(None)
    xlist = []
    ylist = []
    namelist = []
    for key in dictionary:
        xlist.append(dictionary[key].length())
        ylist.append(dictionary[key].gccount())
        namelist.append(dictionary[key].header())

    def onpick(event):
        ind = event.ind
        contigname = namelist[int(ind[0])]   # Print the contig name for 
                                             # the chosen data point. But 
                                             # for cases where two or more
                                             # data points are chosen 
                                             # because they are close 
                                             # together in the plot, it 
                                             # only prints one. But that is
                                             # less of a problem now that
                                             # the picking tolerance is set to
                                             # a very low value (0.5).
        contig.name = contigname
        mark, = ax.plot(xlist[ind[0]], ylist[ind[0]], color = 'r', marker = 'o')
        annotation = ax.text(0.9, 0.05, contig.name, 
                             horizontalalignment='center', 
                             verticalalignment='center', 
                             transform = ax.transAxes)
        fig.canvas.draw()
        annotation.remove()
        ax.lines.remove(mark)
#        print 'Number', ind, namelist[int(ind[0])], np.take(xlist, ind), \
#              np.take(ylist, ind)	# Just a row used for checking that 
                                   	# the right stuff is printed

    # format_coord shows contig.name (if any) in the lower right corner of 
    # the plot, where x and y coordinates are shown. However, as it is now,
    # the cursor needs to be moved slightly before the text is updated with 
    # the new contigname..
    def format_coord(x, y):
        if contig.name:
            return 'x=%.4f, y=%.4f, name: %s'%(x, y, contig.name)
        else:
            return 'x=%.4f, y=%.4f'%(x, y)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.scatter(xlist, ylist, picker=0.5) # Choose a higher picker value for
                                          # higher tolerance when picking 
                                          # points in the plot.
    plt.suptitle('GC - Length', fontsize = 20)
    plt.ylabel('GC content (%)')
    plt.xlabel('Length (nt)')
    fig.canvas.mpl_connect('pick_event', onpick)
    ax.format_coord = format_coord
    plt.show()





# Plot GC content against coverage.
def covgcplot(dictionary):
    contig = conname(None)
    namelist = []
    xlist, xlist1, xlist2, xlist3 = [], [], [], []
    ylist, ylist1, ylist2, ylist3 = [], [], [], []
    for key in dictionary:
        if isNaN(dictionary[key].coverage()) == False:
            xlist.append(dictionary[key].coverage())
            ylist.append(dictionary[key].gccount())
            namelist.append(dictionary[key].header())
            if dictionary[key].length() > 100000:
                xlist3.append(dictionary[key].coverage())
                ylist3.append(dictionary[key].gccount())
            elif 10000 <= dictionary[key].length() <= 100000:
                xlist2.append(dictionary[key].coverage())
                ylist2.append(dictionary[key].gccount())
            elif dictionary[key].length() < 10000:
                xlist1.append(dictionary[key].coverage())
                ylist1.append(dictionary[key].gccount())
            else:
                pass
    def onpick(event):
        ind = event.ind
        contigname = namelist[int(ind[0])]
        contig.name = contigname
        mark, = ax.plot(xlist[ind[0]], ylist[ind[0]], color = 'r', marker = 'o')
        annotation = ax.text(0.9, 0.05, contig.name, 
                             horizontalalignment='center', 
                             verticalalignment='center', 
                             transform = ax.transAxes)
        fig.canvas.draw()
        annotation.remove()
        ax.lines.remove(mark)

    def format_coord(x, y):
        if contig.name:
            return 'x=%.4f, y=%.4f, name: %s'%(x, y, contig.name)
        else:
            return 'x=%.4f, y=%.4f'%(x, y)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.canvas.mpl_connect('pick_event', onpick)
    ax.format_coord = format_coord
    series = plt.scatter(xlist, ylist, marker = 'o', picker = 0.5, 
                         facecolors = 'none', edgecolors = 'none')
    series1 = plt.scatter(xlist1, ylist1, color = 'k', 
                          marker = (5, 2, 0), label="<10k bp")
    series2 = plt.scatter(xlist2, ylist2, facecolors = 'none', 
                          edgecolors = 'b', marker = 'o', 
                          label="10-100k bp")
    series3 = plt.scatter(xlist3, ylist3, facecolors = 'none', 
                          edgecolors = 'r', marker = 'o', 
                          label=">100k bp")
    plt.suptitle('Coverage - GC', fontsize = 20)
    plt.ylabel('GC content (%)')
    plt.xlabel('Coverage')
    leg = plt.legend(scatterpoints = 1)
    leg.get_frame().set_alpha(0.5) # Transparent figure legend so that data
                                   # hiding behind it will still be visible.
    plt.show()





def covlenplot(dictionary):
    contig = conname(None)
    namelist = []
    xlist, xlist1, xlist2, xlist3 = [], [], [], []
    ylist, ylist1, ylist2, ylist3 = [], [], [], []
    for key in dictionary:
        if isNaN(dictionary[key].coverage()) == False:
            xlist.append(dictionary[key].coverage())
            ylist.append(dictionary[key].length())
            namelist.append(dictionary[key].header())
            if dictionary[key].gccount() > 55:
                xlist3.append(dictionary[key].coverage())
                ylist3.append(dictionary[key].length())
            elif 40 <= dictionary[key].gccount() <= 55:
                xlist2.append(dictionary[key].coverage())
                ylist2.append(dictionary[key].length())
            elif dictionary[key].gccount() < 40:
                xlist1.append(dictionary[key].coverage())
                ylist1.append(dictionary[key].length())
            else:
                pass

    def onpick(event):
        ind = event.ind
        contigname = namelist[int(ind[0])]
        contig.name = contigname
        mark, = ax.plot(xlist[ind[0]], ylist[ind[0]], 
                            color = 'c', marker = 'o')
        annotation = ax.text(0.9, 0.05, contig.name, 
                             horizontalalignment='center', 
                             verticalalignment='center', 
                             transform = ax.transAxes)
        fig.canvas.draw()
        annotation.remove()
        ax.lines.remove(mark)

    def format_coord(x, y):
        if contig.name:
            return 'x=%.4f, y=%.4f, name: %s'%(x, y, contig.name)
        else:
            return 'x=%.4f, y=%.4f'%(x, y)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.canvas.mpl_connect('pick_event', onpick)
    ax.format_coord = format_coord
    series = plt.scatter(xlist, ylist, marker = 'o', facecolors = 'none', 
                         edgecolors = 'none', picker = 0.5)
    series1 = plt.scatter(xlist1, ylist1, color = 'k', 
                          marker = 'o', label="<40%")
    series2 = plt.scatter(xlist2, ylist2, color = 'b', 
                          marker = 'o', label="40-55%")
    series3 = plt.scatter(xlist3, ylist3, color = 'r', 
                          marker = 'o', label=">55%")
    plt.suptitle('Length - Coverage', fontsize = 20)
    plt.ylabel('Contig length')
    plt.xlabel('Coverage')
    leg = plt.legend(title = '% GC', scatterpoints = 1)
    leg.get_frame().set_alpha(0.5)
    plt.show()





def covhistogram(dictionary):
    histlist = []
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for key in dictionary:
        if isNaN(dictionary[key].coverage()) == False:
            histlist.append(dictionary[key].coverage())
    plt.hist(histlist, bins=np.logspace(0.1, 7, 200)) # These values can be
                                                      # changed to get another
                                                      # range or bin size.
    ax.set_xscale('log')
    plt.suptitle('Coverage histogram', fontsize = 20)
    plt.xlabel('Coverage')
    plt.ylabel('Frequency')
    plt.show()





def lenhistogram(dictionary):
    histlist = []
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for key in dictionary:
        histlist.append(dictionary[key].length())
    plt.hist(histlist, bins=np.logspace(0.1, 7, 200))
    ax.set_xscale('log')
    plt.suptitle('Length histogram', fontsize = 20)
    plt.xlabel('Length')
    plt.ylabel('Frequency')
    plt.show()





# tests for 'nan'
def isNaN(num):
    return num != num





# Print the outputs that were chosen.
def main():
    dictionary = read_file(args.infile)
    if args.coverage:
        read_covfile(args.coverage, dictionary)
    for key in sorted(dictionary):
        if args.header == True:
            print dictionary[key].header(), '\t',
        if args.length == True:
            print dictionary[key].length(), '\t',
        if args.gccontent == True:
            print dictionary[key].gccount(), '\t',
        if args.ncontent == True:
            print dictionary[key].ncontent(), '\t',
        if args.coverage:
            if isNaN(dictionary[key].coverage()):
                print >> sys.stderr, dictionary[key].coverage(),
            else:
                print dictionary[key].coverage(),
        if len(sys.argv) > 2:    # If no flags are given, 
                                 # no line breaks are printed.
            print	# Just there to introduce a line break.
    if args.lengcplot == True:
        lengcplot(dictionary)
    if args.covgcplot == True:
        try:
            covgcplot(dictionary)
        except:
            sys.stderr.write("ERROR: Correct coverage file not supplied?")
    if args.covlenplot == True:
        covlenplot(dictionary)
    if args.covhistogram == True:
        covhistogram(dictionary)
    if args.lenhistogram == True:
        lenhistogram(dictionary)
    args.infile.close()





if __name__ == "__main__":
    main()
