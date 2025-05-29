#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Python > 3.8

@author: Andy Vierstraete

compare sequences from a MinION run and sort them based on similarity

Lightweight, super fast C/C++ library for sequence alignment using edit (Levenshtein) 
distance:
    https://pypi.org/project/edlib/#description
    python3 -m pip install edlib

"""
from Bio import SeqIO
import multiprocessing
from multiprocessing import Process, Lock, Queue
from threading import Thread
import edlib # faster implementation than Levenshtein for seq longer than 200 bp 
import os
import sys
import random
import time
import datetime
import glob
import pickle
import argparse
import urllib.request
import re
from itertools import zip_longest
import gzip

global tempfile, infile, num_seq, saved_comparelist, comparelist 

version = '2025-05-28'  # version of the script
#==============================================================================
def check_version(version):
    try:   
        link = urllib.request.urlopen('https://raw.githubusercontent.com/avierstr/'
                                      'amplicon_sorter/refs/heads/master/amplicon_sorter.py').read()
        # find the version-date part of the last version on the webpage
        datepart = re.compile(r'(version.*?)(\d{4}-\d{2}-\d{2})(.*version of the script)')
        x = datepart.search(str(link))
        # the 2nd group of the search is the date
        latest_version = x.group(2)
        # compare the date of this version with the version on the webpage
        if version < latest_version:
            version_name = 'amplicon_sorter_' + latest_version + '.py' 
            # download latest version
            urllib.request.urlopen('https://raw.githubusercontent.com/avierstr/'
                                   'amplicon_sorter/refs/heads/master/amplicon_sorter.py')
            urllib.request.urlretrieve('https://raw.githubusercontent.com/avierstr/'
                                       'amplicon_sorter/refs/heads/master/amplicon_sorter.py',
                                       version_name)
            print('\n =====================================================\n'
                  '| NEW VERSION OF AMPLICON_SORTER AVAILABLE            |\n'
                  '| https://github.com/avierstr/amplicon_sorter         |\n'
                  '| Downloaded latest version as:                       |\n' 
                  '|      ' + version_name + '                  |\n'
                  '| Press ctrl-c to exit                                |\n'
                  ' =====================================================\n')
            t = 10
            while t > 0:
                print('Will continue in ' + str(t) + ' seconds...', end='\r')
                time.sleep(1)
                t -= 1
            # to clear previous line completely   
            print('                                                ', end='\r') 
    except:
        pass
#==============================================================================
def get_arguments():
    
    def range_limited_float_type(arg):
        """ Type function for argparse - a float within some predefined bounds """
        try:
            f = float(arg)
        except ValueError:    
            raise argparse.ArgumentTypeError("Must be a floating point number")
        if f < 50 or f > 100:
            raise argparse.ArgumentTypeError("Argument must be > " + str(50.0) 
                                             + " and < " + str(100.0))
        return f

    def range_limited_float_0_200(arg):
        """ Type function for argparse - a float within some predefined bounds """
        try:
            f = float(arg)
        except ValueError:    
            raise argparse.ArgumentTypeError("Must be a floating point number")
        if f < 0 or f > 200:
            raise argparse.ArgumentTypeError("Argument must be > " + str(0.0) 
                                             + " and < " + str(200.0))
        return f

    def valid_file(param): 
        paramlist = [] # make list of files to process
        if os.path.isfile(param): # if input is file
            # check if input file ends on .fastq or .fasta
            base, ext = os.path.splitext(param)
            if ext.lower() not in ('.fasta', '.fastq', '.gz'): 
                raise argparse.ArgumentTypeError('File extension must be .fastq, .fasta or .gz') 
            paramlist.append(param)
        elif os.path.isdir(param): # if input is folder 
            with os.scandir(param) as iterator:
                for file in iterator:
                    if file.name.lower().endswith('.fastq') or file.name.lower().endswith('.fasta') \
                        or file.name.lower().endswith('.gz'):
                        paramlist.append(file.path)
            paramlist.sort()
            if len(paramlist) == 0:
                sys.exit('Can not find files in folder to process.  File extension must be .fastq or .fasta or .gz')
        else:
            sys.exit('Can not find a file or folder to process.  File extension must be .fastq or .fasta or .gz')
        param = paramlist
        return param

    def dir_path(string):
        string = os.path.join(os.getcwd(), string)
        if not os.path.exists(string):
            os.makedirs(string) # create the folder
        return string
            
    parser = argparse.ArgumentParser(description='AmpliconSorter: Sort amplicons\
                                     based on identity and saves them in different\
                                    files including the consensus.' )
    parser.add_argument('-i', '--input', required=True, type = valid_file,
                        help='Input folder or file in fastq, fasta or .gz format')
    parser.add_argument('-min', '--minlength', type = int, required=False, default=300,
                        help='Minimum readlenght to process.  Default=300')
    parser.add_argument('-max', '--maxlength', type = int, required=False, 
                        help='Maximum readlenght to process.  Default=No limit')
    parser.add_argument('-maxr', '--maxreads', type = int, required=False, default=10000,
                        help='Maximum number of reads to process.  Default=10000')
    parser.add_argument('-ar', '--allreads', action = 'store_true',
                        help='Use all reads between length limits to process.\
                            This argument is limited with "-maxreads" to \
                                have a hard limit')
    parser.add_argument('-np', '--nprocesses', type = int, required=False, default=1,
                        help='Number of processors to use. Default=1')
    parser.add_argument('-sg', '--similar_genes', type = range_limited_float_type, 
                        required=False, default=80.0,
                        help='Similarity to sort genes in groups (value between\
                            50 and 100). Default=80.0')
    parser.add_argument('-ssg', '--similar_species_groups', 
                        type = range_limited_float_type, required=False, 
                        help='Similarity to CREATE species groups (value between\
                            50 and 100). Default=estimate value from data')
    parser.add_argument('-ss', '--similar_species', type = range_limited_float_type, 
                        required=False, default=85.0 ,
                        help='Similarity to ADD sequences to a species group\
                            (value between 50 and 100). Default=85.0')
    parser.add_argument('-sc', '--similar_consensus', type = range_limited_float_type, 
                        required=False, default=96.0 ,
                        help='Similarity to COMBINE groups based on the consensus\
                            sequence (value between 50 and 100). Default=96.0')  
    parser.add_argument('-ldc', '--length_diff_consensus', type = range_limited_float_0_200, 
                        required=False, default=8.0 ,
                        help='Length difference (in percent) allowed between consensuses to\
                            COMBINE groups based on the consensus sequence\
                                (value between 0 and 200). Default=8.0')                          
    parser.add_argument('-sfq', '--save_fastq', action = 'store_true',
                        help='Save the results also in fastq files (fastq files\
                            will not contain the consensus sequence)')
    parser.add_argument('-c', '--compressed', action = 'store_true',
                        help='Compress the fasta and fastq files with gz, \
                            not the consensus files')
    parser.add_argument('-ra', '--random', action = 'store_true',
                        help='Takes random reads from the inputfile.')
    parser.add_argument('-a', '--all', action = 'store_true',
                        help='Compare all selected reads with each other.  Only\
                            advised for a small number of reads (< 10000)')
    parser.add_argument('-aln', '--alignment', action = 'store_true',
                        help='Save alignment that is used to create consensus.')  
    parser.add_argument('-amb', '--ambiguous', action = 'store_true',
                        help='Allow search for ambiguous nucleotides.\
                            IUPAC Ambiguity Codes will be used in the consensus')                   
    parser.add_argument('-o', '--outputfolder', type=dir_path, required=False, 
                         help='Save the results in the specified\
                            outputfolder. Default = input folder')
    parser.add_argument('-ho', '--histogram_only', action = 'store_true',
                        help='Only creates a read length histogram.')
    parser.add_argument('-mac', '--macOS', action = 'store_true',
                        help='Option to try when amplicon_sorter crashes on \
                            Mac with M1 processor.')

    args = parser.parse_args()
    return args
#==============================================================================
def save_arguments(): # save all settings in the result.txt file
    outputfolder = args.outputfolder
    with open(os.path.join(outputfolder,'results.txt'), 'a') as rf:
        rf.write('-----------------------------------------------------------\n')
        rf.write('amplicon_sorter version: ' + version + '\n')  
        rf.write('-----------------------------------------------------------\n')
        rf.write('- date and time = ' + datetime.datetime.now().strftime(
            "%B %d, %Y, %I:%M%p") + '\n')    
        rf.write('- input file = ' + os.path.abspath(infolder_file) + '\n')
        # rf.write('- input file = ' + (','.join(args.input)) + '\n')
        rf.write('- output folder = ' + os.path.normpath(args.outputfolder) + '\n')
        rf.write('- minlength = ' + str(args.minlength) + '\n')
        rf.write('- maxlength = ' + str(args.maxlength) + '\n')
        rf.write('- maxreads = ' + str(args.maxreads) + '\n')
        rf.write('- used_reads = ' + '\n')
        rf.write('- allreads = ' + str(args.allreads) + '\n')
        rf.write('- n_processes = ' + str(args.nprocesses) + '\n')
        rf.write('- similar_genes = ' + str(args.similar_genes) + '\n')
        rf.write('- similar_species_groups = ' + str(args.similar_species_groups) 
                 + '\n')
        rf.write('- similar_species = ' + str(args.similar_species) + '\n')
        rf.write('- similar_consensus = ' + str(args.similar_consensus) + '\n')
        rf.write('- length_difference_consensus = ' + str(args.length_diff_consensus) + '\n')
        rf.write('- save_fastq = ' + str(args.save_fastq) + '\n')
        rf.write('- compressed = ' + str(args.compressed) + '\n')
        rf.write('- random = ' + str(args.random) + '\n')
        rf.write('- compare_all = ' + str(args.all) + '\n')
        rf.write('- alignment = ' + str(args.alignment) + '\n')
        rf.write('- ambiguous bases = ' + str(args.ambiguous) + '\n')
        rf.write('- histogram_only = ' + str(args.histogram_only) + '\n')
        rf.write('-----------------------------------------------------------\n')
#==============================================================================  
def distance(X1,X2, mode='NW'):  # calculate the similarity of 2 sequences
    if len(X1) > len(X2): # check which one is longer
        A2 = X1
        A1 = X2
    else:
        A1 = X1
        A2 = X2
    s = edlib.align(A1, A2, task='distance', mode=mode)
    distance = s['editDistance']
    iden = round(1 - distance/len(A2),3)
    return iden
#==============================================================================
def compl_reverse(self):
    inp  = 'ATCGRYKMSW' # translate table for complement
    outp = 'TAGCYRMKSW'
    complement = ''.maketrans(inp, outp)
    R = (self[::-1]).translate(complement)  # complement reverse
    return R
#==============================================================================
def homopolymersort(templist): # sort bases in homopolymer region based on base counts
    templist2 = []
    polymerlist = [] # list for homopolymer regions
    polymerlist.append(templist[0]) # add first base to list
    for x in templist[1:]:
        if x[0][0] == polymerlist[0][0][0]: # if base is the same as the previous one
            polymerlist.append(x)
        else:
            polymerlist.sort(key=lambda x: int(x[0][1]), reverse=True)
            templist2.extend(polymerlist)
            polymerlist = []
            polymerlist.append(x)
    templist2.extend(polymerlist)
    return templist2
#==============================================================================
def degenerate(X1, X2): # make degenerate bases
    if X1[0][0] in ['T', 'C'] and X2[0][0] in ['T', 'C']:
        B = ('Y', X1[0][1] + X2[0][1])
    elif X1[0][0] in ['A', 'G'] and X2[0][0] in ['A', 'G']:
        B = ('R', X1[0][1] + X2[0][1])
    elif X1[0][0] in ['A', 'C'] and X2[0][0] in ['A', 'C']:
        B = ('M', X1[0][1] + X2[0][1])
    elif X1[0][0] in ['G', 'T'] and X2[0][0] in ['G', 'T']:
        B = ('K', X1[0][1] + X2[0][1])
    elif X1[0][0] in ['G', 'C'] and X2[0][0] in ['G', 'C']:
        B = ('S', X1[0][1] + X2[0][1])
    elif X1[0][0] in ['A', 'T'] and X2[0][0] in ['A', 'T']:
        B = ('W', X1[0][1] + X2[0][1])
    elif X1[0][0] in ['Y', 'R', 'M', 'K', 'S', 'W']:
        B = (X1[0][0], X1[0][1] + X2[0][1])
    elif X2[0][0] in ['Y', 'R', 'M', 'K', 'S', 'W']:
        B = (X2[0][0], X1[0][1] + X2[0][1])
    return B
#==============================================================================
def ambiguity(templist, c): # search for ambiguis bases
    templist2 = []
    ambiguitylist = [] # create list for ambiguous bases
    ambiguitylist.append(templist[0]) # add the first base of the list
    for x in templist[1:]:
        if x[0][0] == ambiguitylist[0][0][0]: #if base is the same as the previous one
            if len(x) == 2: # if 2 most abundant bases are arount 100%
                if c*0.35 <= x[0][1] <= c*0.65 and c*0.75 < x[0][1] + x[1][1] < c*1.2:
                    B = degenerate([x[0]], [x[1]])
                    ambiguitylist.append([B])
                else:
                    ambiguitylist.append(x)
            else:
                ambiguitylist.append(x)
            ambiguitylist.sort(key=lambda x: int(x[0][1]), reverse=True)
        else: 
            # ambiguitylist.sort(key=lambda x: int(x[0][1]), reverse=True)
            if c*0.35 <= x[0][1] <= c*0.65 and len(ambiguitylist) < 4: # if the new base is around 50%
                # if new base + last from ambiguity is around 100% 
                if c*0.8 < x[0][1] + ambiguitylist[-1][0][1] < c*1.2:                  
                    B = degenerate(x, ambiguitylist[-1])
                    ambiguitylist[-1] = [B] # replace last base with ambiguis one
                else:
                    templist2.extend(ambiguitylist)
                    ambiguitylist = []
                    # if 2 most abundant bases are arount 100%
                    if len(x) == 2 and c*0.8 < x[0][1] + x[1][1] < c*1.2: 
                        B = degenerate([x[0]], [x[1]])
                        ambiguitylist.append([B])
                    else:
                        ambiguitylist.append(x)
            else:
                templist2.extend(ambiguitylist)
                ambiguitylist = []
                if len(x) == 2:
                    # if 2 most abundant bases are arount 100%
                    if c*0.35 <= x[0][1] <= c*0.65 and c*0.75 < x[0][1] + x[1][1] < c*1.2: 
                        B = degenerate([x[0]], [x[1]])
                        ambiguitylist.append([B])
                    else:
                        ambiguitylist.append(x)
                else:
                    ambiguitylist.append(x)
    templist2.extend(ambiguitylist)
    return templist2
#============================================================================== 
def create_alignment(consensus, readlist):
    # create an alignment out of a list of reads
    alignlist = []
    t = [x for x in consensus] # make list of string
    alignlist.append(t)
    number = re.compile(r'\d+') # the numbers to find in edlib result
    symbol = re.compile(r'\D') # the letters or symbols to find in edlib result
    for x in range(0, len(readlist)):
        q = [b for b in readlist[x]]  
        s = edlib.align(q, t, mode='NW', task='path',
                        additionalEqualities=[("R", "A"), ("R", "G"),
                                              ("Y", "C"), ("Y", "T"),
                                              ("M", "A"), ("M", "C"),
                                              ("K", "G"), ("K", "T"),
                                              ("S", "G"), ("S", "C"),
                                              ("W", "A"), ("W", "T")]) # mode must be NW !!
        scorelist = []
        n = re.findall(number, s['cigar'])
        sy = re.findall(symbol, s['cigar']) 
        for x, y in zip(n, sy):
            scorelist += int(x) * [y]
        insertlist = []
        for i, (x, y) in enumerate(zip_longest(t, scorelist)):
            if y == 'I': # insert gap
                insertlist.append(i)
            if y == 'D': # if delete is needed, insert gap in q
                q.insert(i, '-')
        for i in insertlist:
            for z in alignlist: # if insertion in longest sequence, also 
                                # insert in other aligned sequences
                z.insert(i, '-')
        alignlist.append(q)
    return alignlist
#==============================================================================
def create_consensus(readlist, infile):
    # Make a consensus from a list of reads with the edlib plugin
    aln = args.alignment
    amb = args.ambiguous
    readlist.sort(key=lambda x: len(x), reverse = True)
    readlist2 = readlist.copy() # make copy of list
    consensus = readlist.pop(0) # use longest read as consensus
    for treshold in [0.45, 0.15, 0.5]: 
        alignlist = create_alignment(consensus, readlist)
        # alignlist = realign_front(alignlist)
        # create the consensus
        c = len(alignlist)
        consenlist = []
        templist2 = [] # list with all most abundand bases in alignment
        for y in range(len(alignlist[0])): # for every position in alignment
            tempdict = {}
            for x in alignlist: # for every sequence
                if x[y] != '-':  # for every column in the sequences
                    if x[y] in tempdict:
                        tempdict[x[y]] += 1 
                    else:
                        tempdict[x[y]] = 1
            templist = list(tempdict.items())
            if len(templist) > 0:
                templist.sort(key=lambda x: int(x[1]), reverse=True)
                if templist[0][1] > c*0.10: # if the first base is more than 10% present
                    templist2.append(templist[0:2]) # add the 2 most abundant bases to the list
        templist2 = homopolymersort(templist2) # sort bases in homopolymer region based on base counts
        for x in templist2:
            if x[0][1] > c*treshold:
                consenlist.append(x[0][0])
        consensus = ''.join(consenlist)
        readlist = readlist2 # continue with full list after one cycle
        # random.shuffle(readlist)
    if amb is True: # option to detect ambiguous bases
        templist2 = ambiguity(templist2, c) # search for ambiguis bases
    consenlist = []
    b = 1 # number of bases in homopolymer region
    a = 0 # counting the bases
    # correct for homopolymer region where base count is dropping
    for n, x in enumerate(templist2):
        try:
            if x[0][0] == templist2[n-1][0][0]: # if the previous base is the same
                if x[0][0] in ['A', 'T']:
                    if b >= 4:
                        if x[0][1] > c*0.2:
                            consenlist.append(x[0][0])
                            b += 1
                    else:
                        if x[0][1] > c*treshold:
                            consenlist.append(x[0][0])
                            b += 1
                    a += 1
                elif x[0][0] in ['C', 'G']:
                    if b >= 3:
                        if templist2[n-1][0][1]*0.5 < x[0][1] > c*0.2: # drop can not be more than half of previous base
                            consenlist.append(x[0][0])
                            b += 1
                    else:
                        if x[0][1] > c*treshold:
                            consenlist.append(x[0][0])
                            b += 1
                    a += 1
            else:
                if x[0][1] > c*treshold:
                    consenlist.append(x[0][0])
                    b = 1
                    a += 1
        except IndexError:
            consenlist.append(x[0][0])
    consensus = ''.join(consenlist)
    if aln is True: # option to write alignment to file 
        alignlist = create_alignment(consensus, readlist)
        # outfile = os.path.join(outputfolder, infile.replace('.fasta', '_alignment.fasta').
        #                         replace('.fastq', '_alignment.fasta'))
        # with open(outfile, 'w') as alignm:
        #     # alignm.write(consensus + '\n')
        #     for i, x in enumerate(alignlist):
        #         if i == 0:
        #             i = 'consensus'
        #         seq = ''.join(x)
        #         alignm.write('>' + str(i) + '\n' + seq + '\n')
        #     alignm.write('\n')
    return consensus, alignlist
#==============================================================================
def N50(readlengthlist, bases):  #calculate the N50
    a = sorted(readlengthlist, reverse=True)
    b = int(bases/2)
    N50 = 0
    for x in a:
        N50 += x
        if N50 >= b:
            return x
            break
#==============================================================================
def figure(readlengthlist, total_num_seq):
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib.ticker import AutoMinorLocator
    
    bases = sum(readlengthlist)
    length_N50 = N50(readlengthlist, bases) 
    try:
        max_readlength = sorted(readlengthlist)[-1]
    except IndexError:
        pass
    minlength = args.minlength
    maxlength = args.maxlength
    if maxlength == None:
        maxlength = max_readlength

    plt.figure(1, figsize=[5,5])
    ax = plt.subplot(2,1,1)
    outputfolder = args.outputfolder
    figname = infile.replace('.fastq', '_total_outputfig.pdf').replace('.fasta',
        '_total_outputfig.pdf').replace('.gz', '')
    plt.ylabel('Number of reads')
    plt.title('Read length histogram') 
    plt.text(0.95, 0.55, 'Total yield (Gb): ' + str(round((bases/1000000000),2)) + 
             '\nNumber of reads: ' +
             '{:,}'.format(len(readlengthlist)) + 
             '\n' + str(minlength) + ' < bp < ' + str(maxlength) +
             ': {:,}'.format(total_num_seq) +
             '\nMax readlength: ' + str(round(max_readlength/1000,1)) + ' Kb' \
             '\nN50 = ' + str(round((length_N50/1000),1)) + ' Kb',   
             horizontalalignment='right', transform=plt.gca().transAxes)
    
    plt.hist(readlengthlist, bins='auto', color='green')
    # plot min and max readlength to include
    plt.axvline(minlength, color='red', linewidth=0.8, linestyle='dashed') 
    plt.axvline(maxlength, color='red', linewidth=0.8, linestyle='dashed')

    
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.2f}'.
        format(x/1000))) # divide read length by 1000
    ax.xaxis.set_minor_locator(AutoMinorLocator()) # put subdevisions on x-scale

    ax2 = plt.subplot(2,1,2)
    plt.ylabel('Log Number of reads')
    plt.xlabel('Read length (Kb)')
   
    plt.hist(readlengthlist, log=True, bins='auto', color='green')
    # plot min and max readlength to include
    plt.axvline(minlength, color='red', linewidth=0.8, linestyle='dashed') 
    plt.axvline(maxlength, color='red', linewidth=0.8, linestyle='dashed')
    min_ylim, max_ylim = plt.ylim() # adding text next to min max lines
    plt.text(minlength, 0.7, 'MinLen: {:}'.format(minlength), fontsize=6, 
             rotation=90, horizontalalignment='right')
    plt.text(maxlength, 0.7, 'MaxLen: {:}'.format(maxlength), fontsize=6, 
             rotation=90, horizontalalignment='right')
    
    ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:.2f}'.
        format(x/1000))) # divide read length by 1000
    ax2.xaxis.set_minor_locator(AutoMinorLocator()) # put subdevisions on x-scale
#    plt.legend(loc='center right')
    plt.tight_layout()
    plt.savefig(os.path.join(outputfolder,figname), format='pdf', dpi=300)
    plt.clf() # clear the figure            

    print('Saved "' + figname + '" as a Read Length Histogram.')
#==============================================================================
def read_file(self): # read the inputfile
    # check if it is a gzip file or not, changes the open function
    if infile.endswith('.gz'):
        open_func = gzip.open
    else:  # fasta or fastq
        open_func = open
          
    with open_func(os.path.join(infolder,infile), 'rt') as inf: # check the fileformat
        line = inf.readline()
        if line[0] == '>':
            fileformat = 'fasta'
        elif line[0] == '@':
            fileformat = 'fastq'
            
    minlength = args.minlength
    maxlength = args.maxlength
    maxreads = args.maxreads
    allreads = args.allreads # use all reads
    ran = args.random
    comp_all = args.all # compare all reads with each other
    global comparelist2, num_seq, comparelist
    comparelist = []
    readlengthlist = []
    ti = 0  # total number of reads in the file
    print('Reading ' + self)
    inputfile = open_func(os.path.join(infolder,infile), 'rt')
    for record in SeqIO.parse(inputfile, fileformat):
        ti += 1
        readlengthlist.append(len(record.seq)) # for figure
        # add name, sequence, index number to list
        # index number is needed for resorting to original order later 
        if len(record.seq) >= minlength:
            if maxlength is None:
                comparelist.append([record.id, str(record.seq).upper(), '', '']) 
            else:
                if len(record.seq) <= maxlength:
                    comparelist.append([record.id, str(record.seq).upper(), '', ''])
    total_num_seq = len(comparelist) # total number of reads passed selection
    if total_num_seq < 5: # system hangs with only one sequence, less than 5 has no result
        print('Number of usable sequences is lower than 5, quitting job')
        # sys.exit()
        raise Exception # go to the next file
    for i,[w,x,y,z] in enumerate(comparelist):
        comparelist[i] = [w, x, 'u', i] # replace third pos with 'u' (unique) 
                                        # and 4th pos of each list item with 
                                        # index number
    if allreads == True:
        if total_num_seq <= maxreads:
            maxreads = total_num_seq     
            
    comparelist2 = []
    
    if total_num_seq < 1000:
         comparelist2.append(comparelist)
         sentence = '--> Low number of reads, reading all '
    else:
        if ran == True:
            if comp_all == True:
                if maxreads > total_num_seq:
                    print('Comparing all reads (-a, -all) with each other is not '
                      'possible when selecting more reads than available ('
                      + str(maxreads) + ' selected - ' + str(total_num_seq) + 
                      ' available)')
                    sys.exit()
                else:
                    subcomparelist = random.sample(comparelist, maxreads)
                    comparelist2.append(subcomparelist)
                sentence = '--> Reading random and comparing all '
            else:
                if total_num_seq < 1000:
                    samplesize = total_num_seq
                else:
                    samplesize = 1000
                    
                while maxreads > samplesize:
                    subcomparelist = random.sample(comparelist, samplesize)
                    comparelist2.append(subcomparelist)
                    maxreads = maxreads-samplesize
                else:
                    subcomparelist = random.sample(comparelist, maxreads)
                    comparelist2.append(subcomparelist)
                sentence = '--> Reading random '
        else:
            if comp_all == True:
                if maxreads > total_num_seq:
                    print('Comparing all reads (-a, -all) with each other is not '
                      'possible when selecting more reads than available ('
                      + str(maxreads) + ' selected - ' + str(total_num_seq) + 
                      ' available)')
                    # sys.exit() 
                    raise Exception # go to the next file
                else:
                    comparelist2.append(comparelist[:maxreads])
                    sentence = '--> Reading and comparing all '
            else:
                n = 0
                while maxreads > 1000:
                    subcomparelist = comparelist[n:n+1000]
                    comparelist2.append(subcomparelist) 
                    maxreads = maxreads -1000
                    n += 1000
                else:
                    subcomparelist = comparelist[n:n+maxreads]
                    comparelist2.append(subcomparelist)
                sentence = '--> Reading '
    inputfile.close()
    print(self + ' contains ' + str(ti) + ' reads.')
    
    if args.histogram_only == True: # if only histogram is wanted
        figure(readlengthlist, total_num_seq) # make a read length histogram
        if maxlength is None:
            maxlength = sorted(readlengthlist)[-1] 
        num_seq = sum(len(x) for x in comparelist2) # number of sample reads
        print('--> There are ' + str(total_num_seq) + ' sequences between ' + 
              str(minlength) + ' and ' + str(maxlength) + 'bp')
        # sys.exit()  
        raise Exception # go to the next file 
    else:
        # figure(readlengthlist, total_num_seq) # make a read length histogram 
        num_seq = sum(len(x) for x in comparelist2) # number of sample reads
        if maxlength is None:
            print(sentence + str(num_seq) + ' out of ' + str(total_num_seq) + 
                  ' sequences longer than ' + str(minlength) + 'bp')
        else:
            print(sentence + str(num_seq) + ' out of ' + str(total_num_seq) + 
                  ' sequences between ' + str(minlength) + ' and ' + 
                  str(maxlength) + 'bp')
    return comparelist, comparelist2, num_seq
#==============================================================================
def process_list(self, tempfile): # make files to do comparisons
    global comparelist2, len_todolist
    nprocesses = args.nprocesses
    # compare 1 with 2,3,4,5,...; compare 2 with 3,4,5,... compare 3 with 4,5,...
    chunk = 1000000  # split size of the chunks to feed the multiprocessing queue
    len_todolist = 0
    todoqueue = Queue(maxsize = 2) # max number in queue
    outputfolder = args.outputfolder
    try:  # remove temporary file if exists
        for x in glob.glob(os.path.join(outputfolder, '*.todo')):
            os.remove(x)
        time.sleep(1)
    except FileNotFoundError:
        pass
    
    def queuer(): # function to place items in files for the queue
        global len_todolist, tl
        l = 0  # length of todolist
        k = 0  # number of todofiles
        tl = 0 # total length todolist
        todolist = []
        for d in self: # comparelist2 is a list of lists
            d.sort(key=lambda x: len(x[1])) # sort list based on length seq
            position = 0
            y = 0 #position in first range
            z = 0 #position in 2nd range
            for position in range(position, len(d)-1):
                z = y
                for position2 in range(position+1,len(d)):
                    z +=1
                    A1 = d[position]
                    A2 = d[position2]
                    if len(A1[1])*1.05 < len(A2[1]):
                        pass # don't compare if length of sequences differ to much
                    else:
                        todolist.append([A1,A2])
                        l += 1
                        tl += 1
                        if l == chunk:
                            todofilename = os.path.join(outputfolder, 'file_' + 
                                                        str(k) + '.todo')
                            with open(todofilename, 'wb') as wf:
                                pickle.dump(todolist, wf)
                            todolist = []
                            k += 1
                            l = 0
                y += 1
                for dirpath, dirnames, filenames in os.walk(outputfolder):
                    filenames = [i for i in filenames if i.endswith('.todo')]
                while len(filenames) > 20: # limit numbers of temporary files,
                                           # it can fill a harddisk !
                    time.sleep(10)  
                    for dirpath, dirnames, filenames in os.walk(outputfolder):
                        filenames = [i for i in filenames if i.endswith('.todo')]
        todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
        if tl == 0:
            print('No reads to compare, exiting...')
            with open(os.path.join(outputfolder,'results.txt'), 'a') as rf:
                rf.write('No reads to compare, exiting...')
                tl = -1
        with open(todofilename, 'wb', buffering=0) as wf:
            pickle.dump(todolist, wf)
        todolist = []
        k += 1
        l = 0

        len_todolist = tl
       
        return len_todolist
     
    def feeder(): # function to feed the queue in parts to save memory
        global len_todolist
        while len_todolist == 0:  #start feeding when it is still making todolist
            time.sleep(5)
            for dirpath, dirnames, filenames in os.walk(outputfolder):
                filenames = [i for i in filenames if i.endswith('.todo')]
                filenames.sort(key=lambda x: os.path.getmtime(os.path.join(
                    outputfolder, x))) #sort on modification time
                filenames = filenames[:-1] # don't include last file, possible 
                                           # still being written to disk
                for name in filenames[:]: # only take a few files in memory
                    print('processing: ' + name)
                    with open(os.path.join(outputfolder, name), 'rb') as rf:
                        sublist = pickle.load(rf)
                        todoqueue.put(sublist, block=True) 
                        time.sleep(2)
                    os.remove(os.path.join(outputfolder, name))
        else: # feed all the rest when finished making todolist
            for dirpath, dirnames, filenames in os.walk(outputfolder):
                filenames = [i for i in filenames if i.endswith('.todo')]
                filenames.sort(key=lambda x: os.path.getmtime(os.path.join(
                    outputfolder, x)))
                for name in filenames:
                    print('processing: ' + name)
                    with open(os.path.join(outputfolder, name), 'rb') as rf:
                        sublist = pickle.load(rf)
                        todoqueue.put(sublist, block=True)
                        time.sleep(2)
                    os.remove(os.path.join(outputfolder, name))
            for i in range(nprocesses): # put 'STOP' at the end of the queue 
                                        # for every process
                todoqueue.put("STOP")    
                        
    def consumer(): # function to consume the queue  
        try:
            process = [Process(target=similarity, args=(todoqueue, tempfile,)) 
                        for x in range(nprocesses)]
            for p in process:
                p.start()
            for p in process:
                p.join() 
        except KeyboardInterrupt:
            print("Shutting processes down")
           # Optionally try to gracefully shut down the worker processes here.
            p.terminate()
            p.join()
    
    c = Thread(target = consumer)
    c.start()
    Thread(target = queuer).start()
    time.sleep(5)
    if tl == -1:
        for i in range(nprocesses): # put 'STOP' at the end of the queue 
                                    # for every process
            todoqueue.put("STOP")   
        raise Exception # go to the next file
    Thread(target = feeder).start()
    c.join() # wait until c has finished its work
#==============================================================================
def similarity(todoqueue, tempfile): # process files for similarity 
    global progress, len_todolist
    outputfolder = args.outputfolder
    try:  # remove temporary file if exists
        os.remove(os.path.join(outputfolder, tempfile))
    except FileNotFoundError:
        pass
    similarg = args.similar_genes/100
    templist = []
    MYLOCK = Lock()
    for X in iter(todoqueue.get, 'STOP'):# do stuff until infile.get returns "STOP"
        for A1, A2 in X:
            try:
                b = multiprocessing.current_process()
                iden = distance(A1[1],A2[1])
                if iden >= similarg: 
                    templist.append((str(A1[3]) + ':' + str(A2[3]) + ':' + 
                                     str(iden)))
                elif iden < 0.5:
                    iden = distance(A1[1],compl_reverse(A2[1]))
                    if iden >= similarg:
                        templist.append((str(A1[3]) + ':' + str(A2[3]) + ':' + 
                                         str(iden) + ':' + 'reverse'))
            except KeyboardInterrupt:
                print("Shutting process down")
                b.terminate()
        MYLOCK.acquire()  # save the list every chunk similarities
        with open(os.path.join(outputfolder, tempfile), 'a') as f:
            for c in templist:
                f.write(c +'\n')
        templist =[]
        MYLOCK.release()
#==============================================================================
def SSG(tempfile):  #calculate the N6
    # comparable with the N50 in sequence assembly, here used to estimate ssg value
    outputfolder = args.outputfolder
    print('Estimating the ssg value for this dataset')
    t = 0
    totalsimil = 0
    tempdict = {}
    with open(os.path.join(outputfolder, tempfile), 'r') as tf:
        for line in tf:
            simil = float(line.strip().split(':')[2])
            if simil in tempdict: # if that simil value already exists in dict
                tempdict[simil] += 1
            else:
                tempdict[simil] = 1
            t += 1
            totalsimil += simil # count total value of similarities
    
    templist = list(tempdict.keys())
    templist.sort(reverse=True)
    b = int(totalsimil*0.06) # this value 0.06 looks ok for several testfiles
    N6 = 0
    for x in templist:
        N6 += tempdict[x]*x # simil * number of occurances
        if N6 >= b:
            print('-> Estimated ssg = ' + str(int(x*100)))
            return int(x*100)
            break
#==============================================================================
def finetune(grouplist):
    #--------------------------------------------------------------------------
    def distance_finetune(X1,X2):  # calculate the similarity of 2 sequences
                                    # with the HW mode (begin and end distance
                                    # not important)
        if len(X1) > len(X2): # check which one is longer
            A2 = X1
            A1 = X2
        else:
            A1 = X1
            A2 = X2
        s = edlib.align(A1, A2, task='distance', mode='HW')
        distance = s['editDistance']
        iden = round(1 - distance/len(A2),3) 
        return iden
    #--------------------------------------------------------------------------  
    def reads_direction(readlist): 
        # put all sequences in the same direction (F or R)
        x = readlist[0][1] # first sequence
        for i, y in enumerate(readlist):
            iden = distance(x,y[1])
            idenR = distance(x,compl_reverse(y[1]))
            if iden < idenR:
                readlist[i][1] = compl_reverse(y[1])
        return readlist
    #--------------------------------------------------------------------------
    def check_consensus(consensus, readlist):
        # check if the reads are from one species
        for i, x in enumerate(readlist):  
            iden = distance(consensus,x[1])
            readlist[i][2] = iden
        readlist.sort(key=lambda x: str(x[2]))
        seqlist = [x[1] for x in readlist if x[2] > 0.94]
        if len(seqlist) < 20:
            seqlist = [x[1] for x in readlist[-20:]]
        consensus1, _ = create_consensus(seqlist[-50:], infile)
        iden = distance(consensus1, consensus) 
        return iden, consensus1, readlist
    #--------------------------------------------------------------------------
    addlist = []
    for i, group in enumerate(grouplist):
        group2 = group[:]
        readlist = []
        for n in group:
            if n.isalpha():
                grouplist[i].remove(n)
            else:
                 # get the sequence that matches the number
                 readlist.append([n, comparelist[int(n)][1], '']) 
        
        readlist = reads_direction(readlist) # put reads in same direction
        try:
            readlist2 = random.sample(readlist,100)
        except ValueError:
            readlist2 = random.sample(readlist,len(readlist))
        # compare 100 reads, find a close and distant sequence
        for x in readlist2[0:1]:
            scorelist = []
            for y in readlist2[1:100]:
                iden = distance(x[1], y[1])
                scorelist.append([iden, y[1]])
        scorelist.sort(key=lambda x: x[0])
        consensus1 = scorelist[int(len(scorelist)//1.25)][1] # take seq somewhere at the end
        consensus2 = scorelist[int(len(scorelist)//5)][1] # take seq somewhere at the begin
        
        p = 1
        iden1 = iden2 = 0
        while iden1 < 1 or iden2 < 1:
            iden1, consensus1, readlist = check_consensus(consensus1, readlist)
            iden2, consensus2, readlist = check_consensus(consensus2, readlist)
            iden3 = distance(consensus1, consensus2)
            print('---finetune group ' + str(i) + ' cycle ' + str(p))
            p += 1
            if p == 11:
                break

        if iden3 == 1:
            seqlist = [x[1] for x in readlist if x[2] >= 0.95]
            if len(seqlist) >= 5:
                try:
                    sample = random.sample(seqlist, 150)
                except ValueError:
                    sample = seqlist
                consensus, _ = create_consensus(sample,infile)
                grouplist[i].append(consensus)
                for x in readlist:
                    if x[2] < 0.95: 
                        grouplist[i].remove(x[0])
            else:
                grouplist[i] = []
        else:
            # process group B
            grouplist[i] = [] # empty grouplist
            seqlist = [x[1] for x in readlist if x[2] >= 0.95]
            if len(seqlist) >= 5:
                try:
                    sample = random.sample(seqlist, 150)
                except ValueError:
                    sample = seqlist
                consensus, _ = create_consensus(sample, infile)
                for j, x in enumerate(readlist):
                    if x[2] >= 0.95: 
                        grouplist[i].append(x[0])
                        readlist[j] = []
                grouplist[i].append(consensus)
            if len(grouplist[i]) == 0:
                print('---finetune did not improve it')
                grouplist[i] = group2
            # process group A
            readlist = [x for x in readlist if len(x) > 0]
            if len(readlist) > 5:
                iden1, consensus1, readlist = check_consensus(consensus1, readlist)
                seqlist = [x[1] for x in readlist if x[2] >= 0.95]
                templist = []
                if len(seqlist) >= 5:
                    try:
                        sample = random.sample(seqlist, 150)
                    except ValueError:
                        sample = seqlist
                    consensus, _ = create_consensus(sample, infile)
                    for i, x in enumerate(readlist):
                        if x[2] >= 0.95: 
                            templist.append(x[0])
                    templist.append(consensus)
                    addlist.append(templist)
    if len(addlist) > 0:
        grouplist.extend(addlist)
    return grouplist
#==============================================================================
def update_list(tempfile): # create gene-groups from compared sequences
    outputfolder = args.outputfolder
    if args.similar_species_groups == 'Estimate':
        estimated_ssg = SSG(tempfile) # estimate ssg value
        args.similar_species_groups=estimated_ssg
        with open(os.path.join(outputfolder,'results.txt'), 'r') as rf:
            c = rf.readlines()
            # find position of text
            d = c.index('- similar_species_groups = Estimate\n')
            # replace text
            c[d] = '- similar_species_groups = ' + str(estimated_ssg) + ' (Estimated)\n'
        with open(os.path.join(outputfolder,'results.txt'), 'w') as rf:
            for line in c:
                rf.write(line)
    print('Filtering compared sequences for best hits and create groups')
    global num_seq 
    templist = []
    tempdict = {}
    t = 0 # items done
    try:
        with open(os.path.join(outputfolder, tempfile), 'r') as tf:
            for line in tf:
                e = line.strip().split(':')[:3]
                a = e[1]
                if a in tempdict:
                    tempdict[a].append(e)
                    # sort list based on 2nd number (A2) and score)
                    tempdict[a].sort(key=lambda x: (int(x[1]), float(x[2]))) 
                    for i, j in enumerate(tempdict[a][:-1]):
                        # if second index number is the same for the next item
                        if j[1] == tempdict[a][i+1][1]: 
                            # if iden is lower: keep the best
                            if j[2] < tempdict[a][i+1][2]: 
                                j[0] = ''  # mark to remove
                    # remove those lower values from the list 
                    tempdict[a] = [i for i in tempdict[a] if i[0] != ''] 
                else:
                    tempdict[a] = [e]
                t += 1
                if t % 1000000 == 0: 
                    print(str(t) + ' filtered', end='\r')
        # add all values from dict to a list        
        templist.extend(tempdict.values())
        # make one list out of nested lists
        templist = [x for sublist in templist for x in sublist] 
        # sort list based on score and index number
        templist.sort(key=lambda x: (float(x[2]),int(x[1])),reverse=True) 

    except FileNotFoundError:
        sys.exit()
    
    print('Creating gene groups')
    grouplist = [] 
    q = len(templist)
    r = 0
    # Make groups with sequences with high similarity
    for x in templist:
        r += 1  # give indication how much is done
        if r % 10000 == 0:
            print(str(round(r/q*100, 1)) + '% done', end='\r')
        for s in grouplist: 
            if len({x[0], x[1]}.intersection(s)) > 0:
                s.update({x[0], x[1]})
                break    
        else:
            grouplist.append({x[0], x[1]})
    
    grouplist = merge_groups(grouplist)

    # extra comparison to check of same genes in files
    grouplist = comp_consensus_groups(grouplist) 

    grouped_seq = 0 # number of sequences in grouplist  
    
    for x in grouplist:
        grouped_seq += len(x)  
        
        
    for j, x in enumerate(grouplist):
        outputfile = os.path.join(outputfolder, infile.replace('.fastq', '_').
                                  replace('.fasta', '_').replace('.gz', '') + str(j) + '.group')
        with open(outputfile, 'a') as outputf:
            for y in x:
                outputf.write(str(y) + '\n')
        print('  ' + str(os.path.split(outputfile)[1]) + ' contains ' + 
              str(len(x)) + ' sequences (' + str(round(len(x)*100/num_seq, 2)) + '%)')
    print(str(grouped_seq) + '/' + str(num_seq) + ' sequences assigned in groups (' 
          + str(round(grouped_seq*100/num_seq, 2)) + '%)')

    return grouplist
#==============================================================================
def merge_groups(grouplist):
    a1 = len(grouplist)
    a2 = 0
    if a1 > 1:
        print('--> Number of groups before merge: ' + str(a1)) 
        # remove empty groups who are pointing to an other group from "compare_consensus_groups"
        grouplist = [i for i in grouplist if len(i) > 1]  
        grouplist = [set(i) for i in grouplist]  # make set of the grouplists
        while a1 > a2:
            a1 = len(grouplist) 
            y = 0 #position in first range
            z = 0 #position in 2nd range
            position = 0
            for position in range(position, len(grouplist)-1):
                z = y
                for position2 in range(position+1,len(grouplist)):
                    z +=1
                    A1 = grouplist[position]
                    A2 = grouplist[position2]
                    # check if numbers occur in other groups
                    if len(A1.intersection(A2)) > 0: 
                        grouplist[position] = grouplist[position].union(A2)
                        grouplist[position2].clear() # mark for removal
                y += 1 
            # remove empty subsets
            grouplist = [i for i in grouplist if len(i) > 0]  
            a2 = len(grouplist)
        print('--> Number of groups after merge: ' + str(a2)) 
    
    return grouplist
#==============================================================================
def make_consensus(todoqueue, outputfolder, consensus_tempfile, infile):
    global comparelist
    aln = args.alignment
    compressed = args.compressed
    MYLOCK = Lock()
    templist = []
    for X in iter(todoqueue.get, 'STOP'): # do stuff until infile.get returns "STOP"
        for i, x in X:
            consensuslist = []
            # get the sequence that matches the number
            for y in x:
                consensuslist.append(comparelist[int(y)][1]) 
            consensuslist.sort(key=lambda x: len(x)) #sort list based on length seq
            # get all seq in same direction 
            consensuslist2 = consensus_direction(consensuslist) 
            consensus, alignlist = create_consensus(consensuslist2, infile)  # create consensus sequence
            templist.append([i, consensus])
            if aln is True: # option to write alignment to file 
                if compressed is True: # output files compressed or not
                    open_func = gzip.open
                else:  # fasta or fastq
                    open_func = open
                # try:  # remove temporary alignment files if exists
                #     filename = os.path.join(outputfolder, infile).replace('.group', str(i) + '_alignment.fasta')
                #     for x in glob.glob(os.path.join(outputfolder,filename)):
                #         os.remove(x)      
              
                #     time.sleep(1)
                # except FileNotFoundError:
                #     pass

                if infile != '_':
                    outfile = os.path.join(outputfolder, infile).replace(
                        '.group', '_') + str(i) + '_alignment.fasta'
                    if compressed is True:
                        outfile = outfile + '.gz'
                    with open_func(outfile, 'wt') as alignm:
                        # alignm.write(consensus + '\n')
                        for i, x in enumerate(alignlist):
                            if i == 0:
                                i = 'consensus'
                            seq = ''.join(x)
                            alignm.write('>' + str(i) + '\n' + seq + '\n')
                
    MYLOCK.acquire()
    with open(os.path.join(outputfolder, consensus_tempfile), 'a') as f:
        for c in templist:
            f.write(str(c[0]) + ',' + str(c[1]) +'\n')
        f.flush()
    MYLOCK.release()
#==============================================================================
def iden_consensus(todoqueue, outputfolder, consensus_tempfile, _):
    MYLOCK = Lock()
    templist = []
    for X in iter(todoqueue.get, 'STOP'): # do stuff until infile.get returns "STOP"
        for A1, A2, y, z in X:
            idenlist = []
            iden = distance(A1,A2, mode='HW')
            idenlist.append(iden) # add iden to list
            idenR = distance(A1,compl_reverse(A2), mode='HW')
            idenlist.append(idenR) # add idenR to list
            idenlist.sort(reverse=True) # sort the list
            iden = idenlist[0] # take the biggest value
            if iden >= 0.60: 
                templist.append([y, z, iden])    
    MYLOCK.acquire()
    with open(os.path.join(outputfolder, consensus_tempfile), 'a') as f:
        for c in templist:
            f.write(str(c[0]) + ',' + str(c[1]) + ',' + str(c[2]) +'\n')
        f.flush()
    MYLOCK.release()
#==============================================================================
def do_parallel(outputfolder, nprocesses, consensus_tempfile, make_consensus, 
                  stringx, group_filename):
    for dirpath, dirnames, filenames in os.walk(outputfolder):
        filenames = [i for i in filenames if i.endswith('.todo')]
        filenames.sort(key=lambda x: os.path.getmtime(os.path.join(outputfolder, x)))
        n = 1
        for name in filenames:
            print(stringx + name)
            n += 1
            with open(os.path.join(outputfolder, name), 'rb') as rf:
                todolist = pickle.load(rf)
                if len(todolist) >= nprocesses:
                    chunk = len(todolist)//nprocesses
                else:
                    chunk = 1 #len(todolist)
                todoqueue = Queue()
                b = 0
                e = chunk
                while e < len(todolist):
                    sublist = todolist[b:e]
                    todoqueue.put(sublist)
                    b = e
                    e += chunk
                else:
                    sublist = todolist[b:e]
                    todoqueue.put(sublist)
                time.sleep(2)
            os.remove(os.path.join(outputfolder, name))
            try:
                process = [Process(target=make_consensus, args=(todoqueue, 
                        outputfolder, consensus_tempfile, group_filename,)) for x in range(nprocesses)]
                for p in process:
                    # ask the processes to stop when all files are handled
                    # "STOP" is at the very end of queue
                    todoqueue.put("STOP")
                for p in process:
                    p.start()
                for p in process:
                    p.join()
            except KeyboardInterrupt:
                print("Shutting processes down")
               # Optionally try to gracefully shut down the worker processes here.
                p.terminate()
                p.join()
#==============================================================================
def comp_consensus_groups(grouplist): # compare consensuses with each other
    global comparelist
    comparelist.sort(key=lambda x: x[3]) # sort list based on index number
                                         # must be done for option '-all' sequences !
    outputfolder = args.outputfolder  
    nprocesses = args.nprocesses        
    length_diff_c = args.length_diff_consensus/100 + 1 # length difference allowed for 
                                                # consensus                     
    consensus_tempfile = os.path.join(outputfolder, 'consensus.tmp')
    
    try:  # remove temporary file if exists
        for x in glob.glob(os.path.join(outputfolder, '*.todo')):
            os.remove(x)
        time.sleep(1)
    except FileNotFoundError:
        pass
    try:  # remove temporary file if exists
        os.remove(os.path.join(outputfolder, 'consensus.tmp'))
    except FileNotFoundError:
        pass
                                                                        
    print('-> Merging based on consensus of 50 reads per group') 
    grouplist = [list(i) for i in grouplist]  # make list of the groupsets
    a1 = len(grouplist)
    a2 = 0
    while a1 > a2: # limit number of cycles 
        todolist = []
        a1 = len(grouplist)
        position = 0
        k = 0 # number of todofiles
        for i, x in enumerate(grouplist): # create consensus of 50 seq in each group
            if len(x) > 50: # if number of reads > 50, only take 50 for consensus
                x = random.sample (x, 50)
            todolist.append([i, x])
            if len(todolist) == 500: # save in chuncks to save memory
                todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
                with open(todofilename, 'wb') as wf:
                    pickle.dump(todolist, wf)
                todolist = []
                k += 1
        todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
        with open(todofilename, 'wb') as wf:
            pickle.dump(todolist, wf)
           
        try:  # remove temporary file if exists
            os.remove(os.path.join(outputfolder, 'consensus.tmp'))
        except FileNotFoundError:
            pass    

        stringx = '...making consensuses '
        do_parallel(outputfolder, nprocesses, consensus_tempfile, make_consensus, 
                      stringx, '_')
        
        try:
            with open(consensus_tempfile, 'r') as tf:
                temp = tf.readlines()
            for line in temp:
                i, consensus = line.strip().split(',')
                grouplist[int(i)].append(consensus)  # add consensus sequence to the group 
            os.remove(os.path.join(outputfolder, 'consensus.tmp'))
        except FileNotFoundError:
            pass        
        #---------------------------------------------------------------------  
        todolist = []
        a1 = len(grouplist)
        position = 0
        k = 0 # number of todofiles
        l = 0 # number of comparisons to do
        y = 0 #position in first range
        z = 0 #position in 2nd range
        for position in range(position, len(grouplist)-1):
            z = y
            for position2 in range(position+1,len(grouplist)):
                z +=1
                A1 = grouplist[position][-1]
                A2 = grouplist[position2][-1]
                if len(A1)*length_diff_c < len(A2) or len(A2)*length_diff_c < len(A1):
                    pass  #don't compare if length of sequences differ to much
                else:
                    todolist.append([A1, A2, y, z])
                    l +=1
                    if len(todolist) == 2000000: # save in chuncks to save memory
                        todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
                        with open(todofilename, 'wb') as wf:
                            pickle.dump(todolist, wf)
                        todolist = []
                        k += 1                    
            y += 1          
        todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
        with open(todofilename, 'wb') as wf:
            pickle.dump(todolist, wf)  
        
        stringx = '...comparing consensuses '
        do_parallel(outputfolder, nprocesses, consensus_tempfile, iden_consensus, 
                      stringx, '_')
       
        try:
            with open(consensus_tempfile, 'r') as tf:
                temp = tf.readlines()
                grouplist = [set(x) for x in grouplist]
            for line in temp:
                y, z, iden = line.strip().split(',')
                if length_diff_c > 1.08: # for larger length diff: lot of groups can be merged...
                    if float(iden) >= 0.80: # ... if iden > 0.6
                        while len(grouplist[int(y)]) == 1: # check if groups points to an other group
                            y = int(list(grouplist[int(y)])[0].replace('=', '')) # replace y to the other group
                        grouplist[int(y)].update(grouplist[int(z)]) 
                        grouplist[int(z)] = {'=' + str(y)}
                else:
                    while len(grouplist[int(y)]) == 1: # check if groups points to an other group
                        y = int(list(grouplist[int(y)])[0].replace('=', '')) # replace y to the other group
                    grouplist[int(y)].update(grouplist[int(z)]) 
                    grouplist[int(z)] = {'=' + str(y)}
            os.remove(os.path.join(outputfolder, 'consensus.tmp'))
        except FileNotFoundError:
            pass  

        # # remove empty groups who are pointing to an other group
        # grouplist = [i for i in grouplist if len(i) > 1]  
        
        # merge again based on numbers occuring in different groups 
        grouplist = merge_groups(grouplist) 
    
        for i,x in enumerate(grouplist):  # remove consensus from groups
            grouplist[i] = [y for y in x if y.isdigit()]
        
        a2 = len(grouplist)
    # only keep groups with more than 5 seq
    grouplist = [list(set(i)) for i in grouplist if len(i) > 5] 
    a2 = len(grouplist)
    print('--> Number of groups after removing groups with less than 5 sequences: '
          + str(a2))         
    
    return grouplist
#==============================================================================
def read_indexes(group_filename): # read index numbers from the the group file
    global indexes, comparelist   
    nprocesses = args.nprocesses
    indexes = set()
    outputfolder = args.outputfolder
    consensus_tempfile = os.path.join(outputfolder, 'consensus.tmp')
    similar_species_groups = args.similar_species_groups/100
    try:
        with open(os.path.join(outputfolder, group_filename), 'r') as gf:
            temp = gf.readlines()
            for line in temp:
                indexes.add(line.strip())
        print('reading ' + group_filename + ' containing ' + str(len(indexes)) 
              + ' sequences.')
    except FileNotFoundError:
        print('Can not find ' + group_filename)
        sys.exit()
    if group_filename.endswith('nogroup.group'):
        grouplist = []
    else:
        templist = []
        tempdict = {}
        t = 0 # items done
        try:
            with open(os.path.join(outputfolder, tempfile), 'r') as tf:
                for line in tf:
                    e = line.strip().split(':')
                    a = e[1]
                    if float(e[2]) >= similar_species_groups: 
                        if len({e[0], e[1]}.intersection(indexes)) > 0:
                            if a in tempdict:
                                tempdict[a].append(e)
                                # sort list based on 2nd number (A2) and score)
                                tempdict[a].sort(key=lambda x: (int(x[1]), float(x[2]))) 
                                for i, j in enumerate(tempdict[a][:-1]):
                                    # if second index number is the same for the next item
                                    if j[1] == tempdict[a][i+1][1]: 
                                        # if iden is lower: keep the best
                                        if j[2] < tempdict[a][i+1][2]: 
                                            j[0] = ''  # mark to remove
                                # remove those lower values from the list 
                                tempdict[a] = [i for i in tempdict[a] if i[0] != ''] 
                            else:
                                tempdict[a] = [e]
                            t += 1
                            if t % 1000000 == 0: 
                                print(str(t) + ' processed', end='\r')
            # add all values from dict to a list        
            templist.extend(tempdict.values())
            # make one list out of nested lists
            templist = [x for sublist in templist for x in sublist] 
    
        except FileNotFoundError:
            pass
    
        grouplist = [] 
        # sort list based on score
        templist.sort(key=lambda x: float(x[2]),reverse=True) 
        # Make groups with sequences with high similarity
    
        for x in templist:  
            for s in grouplist: 
                if len({x[1], x[0]}.intersection(s)) > 0:
                    s.update([x[0], x[1]])
                    break    
            else:
                grouplist.append({x[0], x[1]})
    
        grouplist = merge_groups(grouplist)  
        # only keep groups with more than 5 seq
        grouplist = [list(set(i)) for i in grouplist if len(i) > 3] 
        a2 = len(grouplist)
        print('--> Number of groups after removing groups with less than 4 sequences: '
              + str(a2)) 
        
        print('----> Making consensus for each group')
        
        comparelist.sort(key=lambda x: x[3]) # sort list based on index number
                                             # must be done for 'all' sequences !
        try:  # remove temporary file if exists
            for x in glob.glob(os.path.join(outputfolder, '*.todo')):
                os.remove(x)
            time.sleep(1)
        except FileNotFoundError:
            pass
        try:  # remove temporary file if exists
            os.remove(os.path.join(outputfolder, 'consensus.tmp'))
        except FileNotFoundError:
            pass
        
        todolist = []
        k = 0 # number of todofiles
        for i,x in enumerate(grouplist):
            if len(x) > 100: # if number of reads > 100, only take 100 for consensus
                x = random.sample (x, 100)
            todolist.append([i, x])
            if len(todolist) == 500: # save in chuncks to save memory
                todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
                with open(todofilename, 'wb') as wf:
                    pickle.dump(todolist, wf)
                todolist = []
                k += 1
        todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
        with open(todofilename, 'wb') as wf:
            pickle.dump(todolist, wf)

        stringx = '...making consensuses '
        do_parallel(outputfolder, nprocesses, consensus_tempfile, make_consensus, 
                      stringx, group_filename)  
    
        try:
            with open(consensus_tempfile, 'r') as tf:
                temp = tf.readlines()
            for line in temp:
                i, consensus = line.strip().split(',')
                grouplist[int(i)].append(consensus)  # add consensus sequence to the group 
            os.remove(os.path.join(outputfolder, 'consensus.tmp'))
        except FileNotFoundError:
            pass

    return indexes, grouplist
#==============================================================================
def filter_seq(group_filename, grouplist, indexes):
    # filter sequences: put the sequences with high similarity in separate files.
    # Sequences of the same species with the same gene should be in one file.
    if infile.endswith('.gz'): # check if gzip or not
        open_func = gzip.open
    else:  # fasta or fastq
        open_func = open
    with open_func(os.path.join(infolder, infile), 'rt') as inf: # check the fileformat
        line = inf.readline()
        if line[0] == '>':
            fileformat = 'fasta'
        elif line[0] == '@':
            fileformat = 'fastq'
    outputfolder = args.outputfolder    
    compressed = args.compressed
    
    if compressed is True: # output files compressed or not
        open_func = gzip.open
    else:  # fasta or fastq
        open_func = open

    if fileformat == 'fasta': # if the inputfile was fasta, it is not possible to 
        fq = False            # save results in fastq format
    else:
        fq = args.save_fastq # check if it needs to be saved in fastq format
    MYLOCK = Lock()
    print('Writing sequences with high similarity in separate files')
    global comparelist2 
    MYLOCK.acquire()
    if fq == True:
        # index the input fastq file
        try:
            if infile.endswith('.gz'): # is it is a gz file, SeqIO.index can not read it, decompress first
                decompressfile = infile.replace('.gz', '')
                try:
                    record_dict = SeqIO.index(os.path.join(outputfolder, decompressfile), 'fastq')
                except FileNotFoundError:
                    with gzip.open(os.path.join(infolder, infile), 'rt') as zf:
                        with open(os.path.join(outputfolder, decompressfile), 'wt') as of: 
                            d = zf.read(1024)
                            while d:
                                of.write(d)
                                d = zf.read(1024)
                    record_dict = SeqIO.index(os.path.join(outputfolder, decompressfile), 'fastq')
                # os.remove(os.path.join(outputfolder, decompressfile))
            else:
                record_dict = SeqIO.index(os.path.join(infolder, infile), 'fastq')
        except ValueError as e:
            print(e)
    for rec_id, seq, scores, index in comparelist2:
        if str(index) in indexes:
            if scores == 'u':  # sequences that have no similarity with others
                 outputfile = os.path.join(outputfolder, group_filename).replace(
                     '.group', '_') + 'unique.fasta' # unique sequences
                 outputfilefq = os.path.join(outputfolder,group_filename).replace(
                     '.group', '_') + 'unique.fastq' # unique sequences
            else:
                 outputfile = os.path.join(outputfolder,group_filename).replace(
                     '.group', '_') + str(scores) + '.fasta'
                 outputfilefq = os.path.join(outputfolder,group_filename).replace(
                     '.group', '_') + str(scores) + '.fastq'
            if compressed is True:
                outputfile = outputfile + '.gz'
                outputfilefq = outputfilefq + '.gz'
            with open_func(outputfile, 'at') as outputf:
                x = str(seq)
                outputf.write('>' + str(index)  + '\n' + x + '\n')
                if fq == True:
                    with open_func(outputfilefq, 'at') as writer: 
                        SeqIO.write(record_dict[rec_id], writer, 'fastq')
    MYLOCK.release()
    grouped_seq = 0 # number of sequences in grouplist  
    for x in grouplist:
        for y in set(x):
            if y.isdigit(): # list contains consensus sequence, don't count those
                grouped_seq += 1
                
    # consensusfilename = os.path.join(outputfolder, group_filename).replace(
    #     '.group', '_consensussequences.fasta') # group consensusfile
    # total consensusfile
    consensusfile = os.path.join(outputfolder, infile).replace('.fasta', 
    '_consensussequences.fasta').replace('.fastq', '_consensussequences.fasta').replace('.gz', '')
    
    # try:  # remove  file if exists
    #     os.remove(os.path.join(outputfolder, consensusfilename))
    # except FileNotFoundError:
    #     pass 

    MYLOCK.acquire()       
    for j, x in enumerate(grouplist):
        outputfile = os.path.join(outputfolder, group_filename).replace(
            '.group', '_') + str(j) + '.fasta'
        if compressed is True:
            outputfile = outputfile + '.gz'
        consensusname = group_filename.replace('.group', '_') + str(j) 
        
        t = 0
        for y in set(x):
            if y.isalpha(): # if it is a sequence
                with open_func(outputfile, 'at') as outputf:
                    outputf.write('>consensus' + '\n' + y + '\n')
                # with open(consensusfilename, 'a') as outputf:
                #     outputf.write('>consensus_' + str(consensusname) + '(' +
                #                   str(len(x)-1) + ')\n' + y + '\n')
                with open(consensusfile, 'a') as outputf:
                    outputf.write('>consensus_' + str(consensusname) + '(' + 
                                  str(len(x)-1) + ')\n' + y + '\n')
                with open(os.path.join(init_outputfolder, 'consensusfile.fasta'), 'a') as outpf:
                    outpf.write('>consensus_' + str(consensusname) + '(' + 
                                  str(len(x)-1) + ')\n' + y + '\n')
                    
                with open(os.path.join(init_outputfolder, 'results.csv'), 'r') as rc:
                    l = rc.readlines()
                    col = l[0].count(',') # count number of files in first line
                with open(os.path.join(init_outputfolder, 'results.csv'), 'w') as rc:
                    rc.writelines(l)
                    rc.write(consensusname + col*', ' + str(len(x)-1) + '\n')

            elif y.isdigit():
                t += 1 # count number of sequences in group
        with open(os.path.join(outputfolder, 'results.txt'), 'a') as rf:
            try:
                rf.write('--> ' + str(os.path.split(outputfile)[1]) + ' contains ' 
                         + str(t) + ' sequences (' + 
                         str(round(t*100/len(comparelist2), 2)) + '% of total)\n')
            except ZeroDivisionError:
                pass
    with open(os.path.join(outputfolder,'results.txt'), 'a') as rf:
        try:
            if group_filename.endswith('nogroup.group'):
                rf.write(str(len(indexes)) + 
                         ' sequences were not assigned in groups and saved in ' + 
                         str(os.path.split(outputfile)[1]) + ' (' + 
                         str(round(len(indexes)*100/len(comparelist2), 2)) + '% of total)\n')
            else:
                rf.write(str(grouped_seq) + '/' + str(len(indexes)) + 
                         ' sequences assigned in group ' + group_filename + 
                         ' (' + str(round(grouped_seq*100/len(indexes), 2)) + '% of group)\n')
        except ZeroDivisionError:
            pass
    MYLOCK.release()
    
    try:  # remove temporary file if exists
        for x in glob.glob(os.path.join(outputfolder, '*.todo')):
            os.remove(x)
        time.sleep(1)
    except FileNotFoundError:
        pass
    try:  # remove temporary alignment files for unique files if exists
        for dirpath, dirnames, filenames in os.walk(outputfolder):
            for name in filenames:
                if compressed is True:
                    end = '_alignment.fasta.gz'
                    b, a = '_alignment.fasta.gz', '.fasta.gz' # before, after
                else:
                    end = '_alignment.fasta'
                    b, a = '_alignment.fasta', '.fasta'
                if name.endswith(end):
                    partname = name.replace(b, a)
                    if partname not in filenames:
                        os.remove(os.path.join(outputfolder, name))
        time.sleep(1)
    except FileNotFoundError:
        pass
#==============================================================================
def process_consensuslist(indexes, grouplist, group_filename):  
    # compare each sequence with consensus1, consensus2,...
    global len_todolist, comparelist2
    todolist = []
    consensuslist = []
    outputfolder = args.outputfolder
    nprocesses = args.nprocesses
    group_tempfile = os.path.join(outputfolder, group_filename).replace('.group', '.tmp')
    try:  # remove temporary file if exists
        os.remove(group_tempfile)
    except FileNotFoundError:
        pass
    
    try:  # remove temporary file if exists
        for x in glob.glob(os.path.join(outputfolder, '*.todo')):
            os.remove(x)
        time.sleep(1)
    except FileNotFoundError:
        pass   
    
    indexes2 = indexes.copy()  # need a duplicate of indexes to remove items
    for x in grouplist:
        for y in x:
            if y.isdigit(): # list contains consensus sequence, don't check that one
                indexes2.discard(y)  # remove those that are already in a subgroup
    # put a number to each consensus that corresponds to the group            
    for x, y in enumerate(grouplist): 
        consensuslist.append([x, y[-1]])  # add number and consensussequence
    # only keep those from group we are working with
    comparelist4 = [i for i in comparelist2 if str(i[3]) in indexes2] 
    k = 0 # number of todofiles
    l = 0 # number of comparisons to do
    for x in range(0, len(comparelist4)):
        for y in range(0,len(consensuslist)):
            A1 = comparelist4[x]
            A2 = consensuslist[y]
            if len(A1[1])*1.05 < len(A2[1]) or len(A2[1])*1.05 < len(A1[1]):
                pass
            else:
                todolist.append([A1,A2])
                l +=1
                if len(todolist) == 2000000: # save in chuncks to save memory
                    todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
                    with open(todofilename, 'wb') as wf:
                        pickle.dump(todolist, wf)
                    todolist = []
                    k += 1
        if k == 100:
            break
    todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
    with open(todofilename, 'wb') as wf:
        pickle.dump(todolist, wf)

    print(group_filename + '----> ' + str(l) + ' comparisons to calculate')

    try:  # remove temporary file if exists
        os.remove(os.path.join(outputfolder, group_tempfile))
    except FileNotFoundError:
        pass    

    if len(todolist) > 0:
        stringx = '...processing: '
        do_parallel(outputfolder, nprocesses, group_tempfile, similarity_species, 
                      stringx, '_')
#============================================================================== 
def similarity_species(todoqueue, outputfolder, group_tempfile, _):
    # calculate the similarity between sequences
    templist = []
    MYLOCK = Lock()
    for X in iter(todoqueue.get, 'STOP'): # do stuff until infile.get returns "STOP"
        for A1, A2 in X:
            try:
                iden = distance(A1[1],A2[1])
                if iden >= similar - 0.01: # -0.01 to speed up if no seq passed criteria
                    templist.append((str(A1[3]) + ':' + str(A2[0]) + ':' + 
                                     str(iden)))
                elif iden < 0.5 :
                    iden = distance(A1[1],compl_reverse(A2[1]))
                    if iden >= similar -0.01:
                        templist.append((str(A1[3]) + ':' + str(A2[0]) + ':' +
                                         str(iden)))
            except KeyboardInterrupt:
                print("Shutting process down")
    MYLOCK.acquire()
    with open(os.path.join(outputfolder,group_tempfile), 'a') as f:
        for c in templist:
            f.write(c +'\n')
        f.flush()
    MYLOCK.release()
#==============================================================================        
def update_groups(group_filename, grouplist):
    # update groups with sequences with high similarity
    global comparelist2
    nprocesses = args.nprocesses   
    templist = []
    tempdict = {}
    templist2 = []
    min_similar = args.similar_species/100
    outputfolder = args.outputfolder
    consensus_tempfile = os.path.join(outputfolder, 'consensus.tmp')  
    group_tempfile = os.path.join(outputfolder, group_filename).replace('.group', '.tmp')
    try:
        with open(group_tempfile, 'r') as tf:
            for line in tf:
                e = line.strip().split(':')
                a = e[0]
                if a in tempdict:
                    tempdict[a].append(e)
                    # sort list based on 1st number and score
                    tempdict[a].sort(key=lambda x: (int(x[0]), float(x[2])))
                    for i, j in enumerate(tempdict[a][:-1]):
                        # if second index number is the same for the next item
                        if j[0] == tempdict[a][i+1][0]: 
                            # if iden is lower: keep the best
                            if j[2] <= tempdict[a][i+1][2]: 
                                j[0] = ''  # mark to remove
                    # remove those lower values from the list 
                    tempdict[a] = [i for i in tempdict[a] if i[0] != ''] 
                else:
                    tempdict[a] = [e]    
        # add all values from dict to a list        
        templist.extend(tempdict.values())
        # make one list out of nested lists
        templist = [x for sublist in templist for x in sublist] 
    except FileNotFoundError:
        pass
        
    templist2 = [i for i in templist if float(i[2]) >= similar] # only keep those 
    if len(templist2) == 0:
        print(group_filename + '----> no sequences passed criteria')
        if similar <= min_similar:
            try:
                os.remove(group_tempfile)
            except FileNotFoundError:
                pass
    else:
        print(group_filename + '----> ' + str(len(templist2)) + 
              ' sequences added to the groups')
                      
        for x in templist2: 
            grouplist[int(x[1])].append(x[0])
        # remove consensus from groups if new seq are added
        for i,x in enumerate(grouplist):  
            if x[-1].isdigit(): # if the last item is a number, sequence has been added
                grouplist[i] = [y for y in x if y.isdigit()]
     
        print(group_filename + '----> Making consensus for each group')
   
        try:  # remove temporary file if exists
            for x in glob.glob(os.path.join(outputfolder, '*.todo')):
                os.remove(x)
            time.sleep(1)
        except FileNotFoundError:
            pass
        try:  # remove temporary file if exists
            os.remove(os.path.join(outputfolder, 'consensus.tmp'))
        except FileNotFoundError:
            pass
   
        todolist = []
        k = 0 # number of todofiles
        for i, x in enumerate(grouplist):
            if x[-1].isdigit(): # if the last item is a number, sequence has been added
                if len(x) > 200: # if number of reads > 200, only take 200 for consensus
                    x = random.sample(x, 200)
                todolist.append([i, x])
                if len(todolist) == 500: # save in chuncks to save memory
                    todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
                    with open(todofilename, 'wb') as wf:
                        pickle.dump(todolist, wf)
                    todolist = []
                    k += 1
        todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
        with open(todofilename, 'wb') as wf:
            pickle.dump(todolist, wf)

        stringx = '...making consensuses '
        do_parallel(outputfolder, nprocesses, consensus_tempfile, make_consensus, 
                  stringx, group_filename)
        
        try:
            with open(consensus_tempfile, 'r') as tf:
                temp = tf.readlines()
            for line in temp:
                i, consensus = line.strip().split(',')
                grouplist[int(i)].append(consensus)  # add consensus sequence to the group 
            os.remove(os.path.join(outputfolder, 'consensus.tmp'))
        except FileNotFoundError:
            pass         
        
        if similar <= min_similar:
            try:
                os.remove(group_tempfile)
            except FileNotFoundError:
                pass
    
    return comparelist2, grouplist, templist, templist2
#==============================================================================
def consensus_direction(consensuslist): 
    # check if all sequences are in the same direction (F or R)
    consensusset = set()
    for x in consensuslist[0:1]:
        for y in consensuslist[1:]:
            iden = distance(x,y)
            idenR = distance(x,compl_reverse(y))
            if iden > idenR:
                consensusset.update([x, y])
            else:
                consensusset.update([x, compl_reverse(y)])
    consensuslist = list(consensusset)
    return consensuslist
#==============================================================================
def compare_consensus(group_filename, grouplist, length_diff_c):
    # compare consensusses with each other
    similar_consensus = args.similar_consensus/100 
    outputfolder = args.outputfolder
    consensus_tempfile = os.path.join(outputfolder, 'consensus.tmp')
    nprocesses = args.nprocesses  
    try:  # remove temporary file if exists
        for x in glob.glob(os.path.join(outputfolder, '*.todo')):
            os.remove(x)
        time.sleep(1)
    except FileNotFoundError:
        pass
    try:  # remove temporary file if exists
        os.remove(os.path.join(outputfolder, 'consensus.tmp'))
    except FileNotFoundError:
        pass
    
    if len(grouplist) > 1:
        a1 = len(grouplist)
        a2 = 0
        b = 0 
        while a1 > a2 and b < 3: # limit number of cycles
            todolist = []
            a1 = len(grouplist)
            position = 0
            k = 0 # number of todofiles
            l = 0 # number of comparisons to do
            y = 0 #position in first range
            z = 0 #position in 2nd range
            for position in range(position, len(grouplist)-1):
                z = y
                for position2 in range(position+1,len(grouplist)):
                    z +=1
                    A1 = grouplist[position][-1]
                    A2 = grouplist[position2][-1]
                    if len(A1)*length_diff_c < len(A2) or len(A2)*length_diff_c < len(A1):
                        pass  #don't compare if length of sequences differ to much
                    else:
                        todolist.append([A1, A2, y, z])
                        l +=1
                        if len(todolist) == 2000000: # save in chuncks to save memory
                            todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
                            with open(todofilename, 'wb') as wf:
                                pickle.dump(todolist, wf)
                            todolist = []
                            k += 1                    
                y += 1          
            todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
            with open(todofilename, 'wb') as wf:
                pickle.dump(todolist, wf)
               
            stringx = '...comparing consensuses '
            do_parallel(outputfolder, nprocesses, consensus_tempfile, iden_consensus, 
                          stringx, '_')        
            
            try:
                with open(consensus_tempfile, 'r') as tf:
                    temp = tf.readlines()
                    grouplist = [set(x) for x in grouplist]
                for line in temp:
                    y, z, iden = line.strip().split(',')
                    if float(iden) >= similar_consensus:
                        while len(grouplist[int(y)]) == 1: # check if groups points to an other group
                            y = int(list(grouplist[int(y)])[0].replace('=', '')) # replace y to the other group
                        grouplist[int(y)].update(grouplist[int(z)]) 
                        list(set(grouplist[int(y)]))
                        grouplist[int(z)] = {'=' + str(y)}
                os.remove(os.path.join(outputfolder, 'consensus.tmp'))
            except FileNotFoundError:
                pass         
    
            grouplist = merge_groups(grouplist)
        
            for i,x in enumerate(grouplist):  # remove consensus from groups
                grouplist[i] = [y for y in x if y.isdigit()]
                
            try:  # remove temporary file if exists
                for x in glob.glob(os.path.join(outputfolder, '*.todo')):
                    os.remove(x)
                time.sleep(1)
            except FileNotFoundError:
                pass
            try:  # remove temporary file if exists
                os.remove(os.path.join(outputfolder, 'consensus.tmp'))
            except FileNotFoundError:
                pass
            
            todolist = []
            k = 0 # number of todofiles
            for i,x in enumerate(grouplist):
                if len(x) > 200: # if number of reads > 200, only take 200 for consensus
                    x = random.sample (x, 200)
                todolist.append([i, x])
                if len(todolist) == 500: # save in chuncks to save memory
                    todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
                    with open(todofilename, 'wb') as wf:
                        pickle.dump(todolist, wf)
                    todolist = []
                    k += 1
            todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
            with open(todofilename, 'wb') as wf:
                pickle.dump(todolist, wf)

            stringx = '...making consensuses '
            do_parallel(outputfolder, nprocesses, consensus_tempfile, make_consensus, 
                          stringx, group_filename)              
                
            try:
                with open(consensus_tempfile, 'r') as tf:
                    temp = tf.readlines()
                for line in temp:
                    i, consensus = line.strip().split(',')
                    grouplist[int(i)].append(consensus)  # add consensus sequence to the group 
                os.remove(os.path.join(outputfolder, 'consensus.tmp'))
            except FileNotFoundError:
                pass        
            
            a2 = len(grouplist)
            b += 1
        
    return grouplist
#==============================================================================   
def rest_reads(indexes, grouplist, group_filename):
    # compare remaining sequences with the consensussequence of the groups,
    # each time with lower similarity, until stable number of sequences in grouplist
    global similar, comparelist2
    MYLOCK = Lock()
    print(group_filename + 
          '--> Comparing the rest of the sequences with the consensus sequences')
    print(group_filename + '----> similarity = ' + str(round(similar, 2)))
    k = 0
    min_similar = args.similar_species/100
    length_diff_c = args.length_diff_consensus/100 + 1 # length difference allowed for 
                                                # consensus in last step 
    process_consensuslist(indexes, grouplist, group_filename)                
    comparelist2, grouplist, templist, templist2 = update_groups(group_filename, 
                                                                 grouplist)
    if len(templist2) > 200: # first 3 cycles need lot of comparisons, try reduce 
                             # that by merging groups if sequences are added
        # compare consensusses with each other
        grouplist = compare_consensus(group_filename, grouplist, 1.08)
    k += 1
    while similar > min_similar: 
        while k <= 2:
            if len(templist2) > 0:
                process_consensuslist(indexes, grouplist, group_filename)
                comparelist2, grouplist, templist, templist2 = update_groups(
                    group_filename, grouplist)
                if len(templist2) > 0 and k < 2:
                    # compare consensusses with each other
                    if len(grouplist) > 1:
                        grouplist = compare_consensus(group_filename, grouplist, 1.08)  
                k += 1
            else: 
                k = 3
        else:
            k = 0
            if similar in [0.94, 0.88]: # perform 2 finetune cycles
                grouplist = finetune(grouplist)
                grouplist = [x for x in grouplist if len(x) > 0]
                process_consensuslist(indexes, grouplist, group_filename)
            similar = round(similar - 0.01, 2) 
            print(group_filename + '----> similarity = ' + str(similar))
            if len(templist) == 0:
                process_consensuslist(indexes, grouplist, group_filename)
                comparelist2, grouplist, templist, templist2 = update_groups(
                    group_filename, grouplist)
                k += 1
            else: 
                comparelist2, grouplist, templist, templist2 = update_groups(
                    group_filename, grouplist)
                k += 1
    if len(grouplist) > 1:
        # compare consensusses with each other 
        grouplist = compare_consensus(group_filename, grouplist, length_diff_c) 

    MYLOCK.acquire()
    comparelist2.sort(key=lambda x: x[3]) #sort list based on index number
    for o, [rec_id, seq, scores, index] in enumerate(comparelist2):
        for i, n in enumerate(grouplist):
            if str(index) in n:  # it belongs to a group
                comparelist2[o][2] = i      
    MYLOCK.release()
    return grouplist

#==============================================================================
def sort_genes(): # read the input file and sort sequences according to gene groups
    global saved_comparelist, comparelist2
    outputfolder = args.outputfolder
    try:
        filename = infile.replace('.fastq', '_*').replace('.fasta', '_*').replace('.gz', '')
        for x in glob.glob(os.path.join(outputfolder, filename)):
            if not x.endswith(".pdf"): # don't remove readlength histogram
                os.remove(x)
        time.sleep(1)
    except FileNotFoundError:
        pass
    
    if args.histogram_only == True: # if only histogram is wanted
        read_file(infile)
    else:
        read_file(infile)
        process_list(comparelist2, os.path.join(outputfolder, tempfile)) 
        comparelist2 = list(set(([tuple(x) for sublist in comparelist2 for x in 
                                  sublist]))) # make list out of list with sublists
        comparelist2 = [list(x) for x in comparelist2]
        comparelist2.sort(key=lambda x: x[3])
        
        with open(os.path.join(outputfolder, saved_comparelist), 'wb') as wf:
            pickle.dump(comparelist, wf)
            pickle.dump(comparelist2, wf)
    # write number of used reads in result file
    with open(os.path.join(outputfolder,'results.txt'), 'r') as rf:
        c = rf.readlines()
        # find position of text
        d = c.index('- used_reads = ' + '\n')
        # replace text
        c[d] = '- used_reads = ' + str(len(comparelist2)) + '\n'
    with open(os.path.join(outputfolder,'results.txt'), 'w') as rf:
        for line in c:
            rf.write(line)
    with open(os.path.join(init_outputfolder, 'results.csv'), 'r') as rc:
        l = rc.readlines()
        l[0] = l[0].strip() + ', ' + filename.replace('_*', '') + '\n' # add name to first line
        col = l[0].count(',')
    with open(os.path.join(init_outputfolder, 'results.csv'), 'w') as rc:
        rc.writelines(l)
        rc.write('Total' + col*', ' + str(len(comparelist)) + '\n')
#==============================================================================
def sort_groups(): # read the gene groups and sort sequences to species level
    global comparelist, comparelist2, resultlist, num_seq, saved_comparelist
    outputfolder = args.outputfolder
    try:
        with open(os.path.join(outputfolder, saved_comparelist), 'rb') as rf:
            comparelist = pickle.load(rf)
            comparelist2 = pickle.load(rf)
        num_seq = len(comparelist2)
    except FileNotFoundError:
        pass 
    
    consensusfile = infile.replace('.fasta', '_consensussequences.fasta').replace(
        '.fastq', '_consensussequences.fasta').replace('.gz', '') # total consensusfile
    try:  # remove  file if exists
        os.remove(os.path.join(outputfolder, consensusfile))
    except FileNotFoundError:
        pass   
    
    grouplist = update_list(os.path.join(outputfolder, tempfile))
    #------------------------------------
    # find numbers of sequences that are not in groups and save in separate file
    groupedseq = [int(x) for sublist in grouplist for x in sublist] # make flat list of numbers
    comparedseq = [x[3] for x in comparelist2]
    unique = set(comparedseq).difference(set(groupedseq))
    if len(unique) > 0:
        outputfile = os.path.join(outputfolder, infile.replace('.fastq', '_').
                                  replace('.fasta', '_').replace('.gz', '') + 'nogroup.group')
        with open(outputfile, 'a') as outputf:
            for y in unique:
                outputf.write(str(y) + '\n')
    #------------------------------------
    for dirpath, dirnames, filenames in os.walk(outputfolder):
        for name in sorted(filenames):
            if name.endswith('.group'):
                sort(name)
    # print info of groups at the end
    with open(os.path.join(outputfolder,'results.txt'), 'r') as rf:
        results = rf.readlines()
        for line in results:
            print(line.strip())
#==============================================================================
def sort(group_filename): 
    global similar, comparelist
    outputfolder = args.outputfolder
           
    try:  # remove temporary file if exists
        filename = os.path.join(outputfolder, group_filename).replace('.group', '_*')
        for x in glob.glob(os.path.join(outputfolder,filename)):
            os.remove(x)      
  
        time.sleep(1)
    except FileNotFoundError:
        pass

    indexes, grouplist = read_indexes(group_filename)
    if len(grouplist) == 0: # if there are no groups, no need to compare rest of sequences
        print('--> Nothing to do')
        filter_seq(group_filename, grouplist, indexes)
        os.remove(os.path.join(outputfolder, group_filename))
    else:
        similar = 0.95
        grouplist = rest_reads(indexes, grouplist, group_filename)
        filter_seq(group_filename, grouplist, indexes)
        os.remove(os.path.join(outputfolder, group_filename))        
#==============================================================================    
if __name__ == '__main__':
    try:
        args = get_arguments()
        if args.macOS is True:
            multiprocessing.set_start_method('fork')
        check_version(version)
        if args.similar_species_groups is None: # needed if multiple files are done
            ssg = 'Estimate'                    # otherwise it uses the same value for all
        else:
            ssg = ''
        init_outputfolder = args.outputfolder
        infolder_file_list = args.input
        for infolder_file in infolder_file_list:
            try:
                if ssg == 'Estimate':
                    args.similar_species_groups = 'Estimate'
                infolder, infile = os.path.split(os.path.realpath(infolder_file))
                tempfile = infile.replace('.fastq', '_compare.tmp').replace('.fasta', 
                                                        '_compare.tmp').replace('.gz', '')
                saved_comparelist = infile.replace('.fastq', '_comparelist.pickle').replace(
                    '.fasta', '_comparelist.pickle').replace('.gz', '')
                outfolder, ext = os.path.splitext(infile.replace('.gz', ''))
                if init_outputfolder:
                    if outfolder not in init_outputfolder:
                        try:
                            args.outputfolder = os.path.join(init_outputfolder, outfolder)
                            os.makedirs(args.outputfolder)
                        except FileExistsError:
                            pass
                else:
                    init_outputfolder = infolder
                    try:
                        args.outputfolder = os.path.join(infolder, outfolder)
                        os.makedirs(args.outputfolder)
                    except FileExistsError:
                        pass
                save_arguments() # write all settings in the results.txt file
                try: # initialise csv result file with one line
                    with open(os.path.join(init_outputfolder, 'results.csv'), 'r') as rc:
                        _ = rc.readline()
                except FileNotFoundError:
                    with open(os.path.join(init_outputfolder, 'results.csv'), 'w') as rc:
                        rc.write(' \n')
                sort_genes()
                sort_groups()
                outputfolder = args.outputfolder
                os.remove(os.path.join(outputfolder, tempfile))
                os.remove(os.path.join(outputfolder, saved_comparelist))
                if infile.endswith('.gz') and args.save_fastq is True: # remove decompressed file from indexing
                    decompressfile = infile.replace('.gz', '')
                    os.remove(os.path.join(outputfolder, decompressfile))
            except Exception: 
                continue
    except KeyboardInterrupt:
        sys.exit()