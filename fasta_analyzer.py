#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
Copyright (C) 2014 Sandra Karlsten. sandra.karlsten@gmail.com

Citation: If you use this version of the program, please cite;
Sandra Karlsten (2014) 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
"""

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





LENGTH_LARGE = 100000
LENGTH_SMALL = 10000
GC_LARGE = 55
GC_SMALL = 40
mark = None
annotation = None





class Fasta(object):
    """Each sequence becomes one object of the class Fasta. 
    """
    def __init__(self, name, seq):
        self.name = name[1:].rstrip() # The name equals the contig header minus the '>'.
        self.seq = seq.lower()
        self.cov = float('nan')	# For missing coverage values, 
                                # covgcplot will not plot.

    def header(self):
        return self.name

    def sequence(self):
        return self.seq

    def length(self):
        return len(self.seq) - self.seq.count('\n')

    # Calculates the number of G and C relative to the 
    # total number of G, C, A and T.
    def gccount(self):
        self.gc = self.seq.count("g") + self.seq.count("c")
        self.total = self.seq.count("g") + self.seq.count("c") + \
                     self.seq.count("a") + self.seq.count("t")
        if self.total == 0:
            self.content = float('nan')
        else:
            self.content = round((float(self.gc) / self.total) * 100, 1)
        return self.content

    def getcoverage(self):
        return self.cov

    def setcoverage(self, cov):
        self.cov = cov

    def ncontent(self):
        return self.seq.count("n")





class conname(object):
    """Store the contig name so that it is easily available to the onpick and 
    format_coord functions."""
    def __init__(self, name):
        self.name = name





def read_covfile(infile, dictionary):
    """Read coverage file and add coverage for each contig to the dictionary.
    The arguments needed are a fasta file and a dictionary where the item is
    an object of the class Fasta and the key is the name of the contig.
    """
    infile.seek(0)
    i = 0
    for line in infile:
        i = i + 1
        name, cov = line.split('\t')
        if name in dictionary.keys():
            try:
                dictionary[name].setcoverage(float(cov))
            except:
                pass
        else:
            sys.stderr.write("ERROR: Contig names did not match, line %r in coverage file. \n" % i)
    return dictionary





def read_file(infile):
    """Read fasta file. The argument needed is a fasta file.
    """
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





def lengcplot(dictionary):
    """Plot GC content against length. The argument needed is a dictionary 
    where the item is an object of the class Fasta and the key is the name 
    of the contig.
    """
    contig = conname(None)
    xlist = []
    ylist = []
    namelist = []
    for key in dictionary:
        xlist.append(dictionary[key].length())
        ylist.append(dictionary[key].gccount())
        namelist.append(dictionary[key].header())

    def onpick(event):
        global annotation
        global mark
        if mark != None: 
            try:
                ax.lines.remove(mark)
                annotation.remove()
            except:
                pass
        """ Print the contig name for the chosen data point. But for cases where
        two or more data points are chosen because they are close together in 
        the plot, it only prints one. But that is less of a problem now that 
        the picking tolerance is set to a very low value (0.5).
        """
        ind = event.ind
        contigname = namelist[int(ind[0])]
        contig.name = contigname
        contig_info = "%s \n Length: %r \n GC: %r \n Coverage: %r \n N: %r" \
                      %(contig.name, dictionary[contig.name].length(), 
                      dictionary[contig.name].gccount(), 
                      dictionary[contig.name].getcoverage(), 
                      dictionary[contig.name].ncontent())
        mark, = ax.plot(xlist[ind[0]], ylist[ind[0]], color = 'r', marker = 'o')
        annotation = ax.text(0.8, 0.15, contig_info, 
                             horizontalalignment='center', 
                             verticalalignment='center', 
                             transform = ax.transAxes)
        def on_key(event):
            global annotation
            global mark
            if event.key == 'control':
                try:
                    annotation.remove()
                    ax.lines.remove(mark)
                    fig.canvas.draw()
                except: 
                    pass
            else:
                pass

        cid = fig.canvas.mpl_connect('key_press_event', on_key)
        fig.canvas.draw()

    def format_coord(x, y):
        """ format_coord shows contig.name (if any) in the lower right corner
        of the plot, where x and y coordinates are shown. However, as it is now,
        the cursor needs to be moved slightly before the text is updated with 
        the new contigname..
        """
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
    return xlist, ylist, namelist





def covgcplot(dictionary):
    """Plot GC content against coverage. The argument needed is a dictionary 
    where the item is an object of the class Fasta and the key is the name 
    of the contig.
    """
    contig = conname(None)
    namelist = []
    xlist, xlist1, xlist2, xlist3 = [], [], [], []
    ylist, ylist1, ylist2, ylist3 = [], [], [], []
    for key in dictionary:
        if isNaN(dictionary[key].getcoverage()) == False:
            xlist.append(dictionary[key].getcoverage())
            ylist.append(dictionary[key].gccount())
            namelist.append(dictionary[key].header())
            if dictionary[key].length() > LENGTH_LARGE:
                xlist3.append(dictionary[key].getcoverage())
                ylist3.append(dictionary[key].gccount())
            elif LENGTH_SMALL <= dictionary[key].length() <= LENGTH_LARGE:
                xlist2.append(dictionary[key].getcoverage())
                ylist2.append(dictionary[key].gccount())
            elif dictionary[key].length() < LENGTH_SMALL:
                xlist1.append(dictionary[key].getcoverage())
                ylist1.append(dictionary[key].gccount())
            else:
                pass
    def onpick(event):
        global annotation
        global mark
        if mark != None: 
            try:
                ax.lines.remove(mark)
                annotation.remove()
            except:
                pass
        ind = event.ind
        contigname = namelist[int(ind[0])]
        contig.name = contigname
        contig_info = "%s \n Length: %r \n GC: %r \n Coverage: %r \n N: %r" \
                      %(contig.name, dictionary[contig.name].length(), 
                      dictionary[contig.name].gccount(), 
                      dictionary[contig.name].getcoverage(), 
                      dictionary[contig.name].ncontent())
        mark, = ax.plot(xlist[ind[0]], ylist[ind[0]], color = 'r', marker = 'o')
        annotation = ax.text(0.8, 0.15, contig_info, 
                             horizontalalignment='center', 
                             verticalalignment='center', 
                             transform = ax.transAxes)
        def on_key(event):
            global annotation
            global mark
            if event.key == 'control':
                try:
                    annotation.remove()
                    ax.lines.remove(mark)
                    fig.canvas.draw()
                except: 
                    pass
            else:
                pass

        cid = fig.canvas.mpl_connect('key_press_event', on_key)
        fig.canvas.draw()

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
                          marker = (5, 2, 0), label="< %s bp"%(LENGTH_SMALL))
    series2 = plt.scatter(xlist2, ylist2, facecolors = 'none', 
                          edgecolors = 'b', marker = 'o', 
                          label="%s - %s bp" %(LENGTH_SMALL, LENGTH_LARGE))
    series3 = plt.scatter(xlist3, ylist3, facecolors = 'none', 
                          edgecolors = 'r', marker = 'o', 
                          label="> %s bp" %(LENGTH_LARGE))
    plt.suptitle('Coverage - GC', fontsize = 20)
    plt.ylabel('GC content (%)')
    plt.xlabel('Coverage')
    leg = plt.legend(scatterpoints = 1)
    leg.get_frame().set_alpha(0.5) # Transparent figure legend so that data
                                   # hiding behind it will still be visible.
    plt.show()
    return xlist, xlist1, xlist2, xlist3, ylist, ylist1, ylist2, ylist3,\
    namelist





def covlenplot(dictionary):
    """Plot coverage against length. The argument needed is a dictionary 
    where the item is an object of the class Fasta and the key is the 
    name of the contig.
    """
    contig = conname(None)
    namelist = []
    xlist, xlist1, xlist2, xlist3 = [], [], [], []
    ylist, ylist1, ylist2, ylist3 = [], [], [], []
    for key in dictionary:
        if isNaN(dictionary[key].getcoverage()) == False:
            xlist.append(dictionary[key].getcoverage())
            ylist.append(dictionary[key].length())
            namelist.append(dictionary[key].header())
            if dictionary[key].gccount() > GC_LARGE:
                xlist3.append(dictionary[key].getcoverage())
                ylist3.append(dictionary[key].length())
            elif GC_SMALL <= dictionary[key].gccount() <= GC_LARGE:
                xlist2.append(dictionary[key].getcoverage())
                ylist2.append(dictionary[key].length())
            elif dictionary[key].gccount() < GC_SMALL:
                xlist1.append(dictionary[key].getcoverage())
                ylist1.append(dictionary[key].length())
            else:
                pass

    def onpick(event):
        global annotation
        global mark
        if mark != None: 
            try:
                ax.lines.remove(mark)
                annotation.remove()
            except:
                pass
        ind = event.ind
        contigname = namelist[int(ind[0])]
        contig.name = contigname
        contig_info = "%s \n Length: %r \n GC: %r \n Coverage: %r \n N: %r" \
                      %(contig.name, dictionary[contig.name].length(), 
                      dictionary[contig.name].gccount(), 
                      dictionary[contig.name].getcoverage(), 
                      dictionary[contig.name].ncontent())
        mark, = ax.plot(xlist[ind[0]], ylist[ind[0]], 
                            color = 'c', marker = 'o')
        annotation = ax.text(0.8, 0.15, contig_info, 
                             horizontalalignment='center', 
                             verticalalignment='center', 
                             transform = ax.transAxes)
        def on_key(event):
            global annotation
            global mark
            if event.key == 'control':
                try:
                    annotation.remove()
                    ax.lines.remove(mark)
                    fig.canvas.draw()
                except: 
                    pass
            else:
                pass

        cid = fig.canvas.mpl_connect('key_press_event', on_key)
        fig.canvas.draw()

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
                          marker = 'o', label="< %r %%" %(GC_SMALL))
    series2 = plt.scatter(xlist2, ylist2, color = 'b', 
                          marker = 'o', label=" %r - %r %%" %(GC_SMALL, GC_LARGE))
    series3 = plt.scatter(xlist3, ylist3, color = 'r', 
                          marker = 'o', label="> %r %%" %(GC_LARGE))
    plt.suptitle('Length - Coverage', fontsize = 20)
    plt.ylabel('Contig length')
    plt.xlabel('Coverage')
    leg = plt.legend(title = '% GC', scatterpoints = 1)
    leg.get_frame().set_alpha(0.5)
    plt.show()
    return xlist, xlist1, xlist2, xlist3, ylist, ylist1, ylist2, ylist3,\
    namelist





def covhistogram(dictionary):
    """Create a histogram over coverage. The argument needed is a dictionary 
    where the item is an object of the class Fasta and the key is the name 
    of the contig.
    """
    histlist = []
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for key in dictionary:
        if isNaN(dictionary[key].getcoverage()) == False:
            histlist.append(dictionary[key].getcoverage())
    plt.hist(histlist, bins=np.logspace(0.1, 7, 200)) # These values can be
                                                      # changed to get another
                                                      # range or bin size.
    ax.set_xscale('log')
    plt.suptitle('Coverage histogram', fontsize = 20)
    plt.xlabel('Coverage')
    plt.ylabel('Frequency')
    plt.show()
    return histlist





def lenhistogram(dictionary):
    """Create a histogram over length. The argument needed is a dictionary 
    where the item is an object of the class Fasta and the key is the name 
    of the contig.
    """
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
    return histlist





def isNaN(num):
    """tests for 'nan'. The argument needed is the number that is to be tested.
    """
    return num != num





def run(args):
    """ Print the outputs that were chosen. The argument needed is the 
    parser.parse_args() that determines which flags and infiles the 
    script can use.
    """
    dictionary = read_file(args.infile)
    if args.coverage:
        read_covfile(args.coverage, dictionary)
    if args.allflags: # If "-all"-flag is given, set flags to True.
        args.header = True
        args.length = True
        args.gccontent = True
        args.ncontent = True
        args.lengcplot = True
        args.lenhistogram = True
        if args.coverage:
            args.covgcplot = True
            args.covlenplot = True
            args.covhistogram = True
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
            if isNaN(dictionary[key].getcoverage()):
                print >> sys.stderr, dictionary[key].getcoverage(),
            else:
                print dictionary[key].getcoverage(),
        if len(sys.argv) > 2:    # If no flags are given, 
                                 # no line breaks are printed.
            print	# Just there to introduce a line break.
    if args.lengcplot == True:
        lengcplot(dictionary)
    if args.covgcplot == True:
        if args.coverage:
            try:
                covgcplot(dictionary)
            except:
                sys.stderr.write("ERROR: Correct coverage file not supplied?")
        else:
            sys.stderr.write("ERROR: Can't run function 'covgcplot'. \
                              No coverage file supplied.\n")
    if args.covlenplot == True:
        if args.coverage:
            covlenplot(dictionary)
        else:
            sys.stderr.write("ERROR: Can't run function 'covlenplot'. \
                              No coverage file supplied.\n")
    if args.covhistogram == True:
        if args.coverage:
            covhistogram(dictionary)
        else:
            sys.stderr.write("ERROR: Can't run function 'covhistogram'. \
                              No coverage file supplied.\n")
    if args.lenhistogram == True:
        lenhistogram(dictionary)
    args.infile.close()
    if args.coverage:
        args.coverage.close()





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 
                                 """Analyze data from a fasta file and the 
                                    corresponding coverage file by calculating
                                    length, GC content and N content and by 
                                    creating scatter plots and histograms. 
                                    """)

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

    parser.add_argument("-all", "--allflags",
                    help = "Shortcut for using all flags.", # Plot functions 
                    action = "store_true")                  # that plot 
                                                            # coverage will 
                                                            # not be used if 
                                                            # no coverage file
                                                            # is provided.

    args = parser.parse_args()
    run(args)
