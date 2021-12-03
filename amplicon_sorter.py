#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Python 3.8

@author: Andy Vierstraete

compare sequences from a MinION run and sort them based on similarity

For c implementation of Levenshtein:
    https://pypi.org/project/python-Levenshtein/
    sudo apt-get install python3-dev
    (needed: python3-setuptools, python3-pip, python3-wheel)
    python3 -m pip install python-Levenshtein
Lightweight, super fast C/C++ library for sequence alignment using edit (Levenshtein) 
distance:
    https://pypi.org/project/edlib/#description
    python3 -m pip install edlib

"""
from Bio import SeqIO
import multiprocessing
from multiprocessing import Process, Lock, Queue
from threading import Thread
# from Levenshtein import distance as l_distance
from Levenshtein import median as l_median
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

global tempfile, infile, num_seq, saved_comparelist, comparelist

version = '2021-12-01'  # version of the script

#==============================================================================
def check_version(version):
    try:   
        link = urllib.request.urlopen('https://github.com/avierstr/amplicon_sorter'
                                      '/blob/master/amplicon_sorter.py').read()
        # find the version-date part of the last version on the webpage
        datepart = re.compile(r'(version</span>.*?)(\d{4}-\d{2}-\d{2})(.*version of the script)')
        x = datepart.search(str(link))
        # the 2nd group of the search is the date
        latest_version = x.group(2)
        # compare the date of this version with the version on the webpage
        if version < latest_version:
            version_name = 'amplicon_sorter_' + latest_version + '.py' 
            # download latest version
            urllib.request.urlopen('https://raw.githubusercontent.com/avierstr/'
                                   'amplicon_sorter/master/amplicon_sorter.py')
            urllib.request.urlretrieve('https://raw.githubusercontent.com/avierstr'
                                       '/amplicon_sorter/master/amplicon_sorter.py',
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

    def valid_file(param):  # check if input file ends on .fastq or .fasta
        base, ext = os.path.splitext(param)
        if ext.lower() not in ('.fasta', '.fastq'): 
            raise argparse.ArgumentTypeError('File extension must be .fastq or .fasta') 
        return param

    def dir_path(string):
        if os.path.exists(string):
            string = os.path.join(os.getcwd(), string)
            return string
        else:
            string = os.path.join(os.getcwd(), string)
            os.makedirs(string) # create the folder
            return string
            
    parser = argparse.ArgumentParser(description='AmpliconSorter: Sort amplicons\
                                     based on identity and saves them in different\
                                    files including the consensus.' )
    parser.add_argument('-i', '--input', required=True, type = valid_file,
                        help='Input file in fastq or fasta format')
    parser.add_argument('-min', '--minlength', type = int, required=False, default=300,
                        help='Minimum readlenght to process.  Default=300')
    parser.add_argument('-max', '--maxlength', type = int, required=False, 
                        help='Maximum readlenght to process.  Default=No limit')
    parser.add_argument('-maxr', '--maxreads', type = int, required=False, default=10000,
                        help='Maximum number of reads to process.  Default=10000')
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
    parser.add_argument('-sfq', '--save_fastq', action = 'store_true',
                        help='Save the results also in fastq files (fastq files\
                            will not contain the consensus sequence)')
    parser.add_argument('-ra', '--random', action = 'store_true',
                        help='Takes random reads from the inputfile.')
    parser.add_argument('-a', '--all', action = 'store_true',
                        help='Compare all selected reads with each other.  Only\
                            advised for a small number of reads (< 10000)')
    parser.add_argument('-o', '--outputfolder', type=dir_path, required=False, 
                        default='./', help='Save the results in the specified\
                            outputfolder. Default = current working directory')
    parser.add_argument('-ho', '--histogram_only', action = 'store_true',
                        help='Only creates a read length histogram.')
    parser.add_argument('-so', '--species_only', action = 'store_true',
                        help='Only create species groups and sort to species level.')

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
        rf.write('- input file = ' + args.input + '\n')
        rf.write('- output folder = ' + args.outputfolder + '\n')
        rf.write('- minlength = ' + str(args.minlength) + '\n')
        rf.write('- maxlength = ' + str(args.maxlength) + '\n')
        rf.write('- maxreads = ' + str(args.maxreads) + '\n')
        rf.write('- n_processes = ' + str(args.nprocesses) + '\n')
        rf.write('- similar_genes = ' + str(args.similar_genes) + '\n')
        rf.write('- similar_species_groups = ' + str(args.similar_species_groups) 
                 + '\n')
        rf.write('- similar_species = ' + str(args.similar_species) + '\n')
        rf.write('- similar_consensus = ' + str(args.similar_consensus) + '\n')
        rf.write('- save_fastq = ' + str(args.save_fastq) + '\n')
        rf.write('- random = ' + str(args.random) + '\n')
        rf.write('- all_reads = ' + str(args.all) + '\n')
        rf.write('- histogram_only = ' + str(args.histogram_only) + '\n')
        rf.write('- species_only = ' + str(args.species_only) + '\n')
        rf.write('-----------------------------------------------------------\n')
#==============================================================================  
def distance(X1,X2):  # calculate the similarity of 2 sequences
    if len(X1) > len(X2): # check which one is longer
        A2 = X1
        A1 = X2
    else:
        A1 = X1
        A2 = X2
    # if A1 is much shorter than A2, it influences the similarity
    if len(A1)*1.05 < len(A2): 
        idenprev = 0
        templist = []
        for a in range(0, len(A2) - len(A1),2): 
            e = int(len(A1)*1.015) # leave room for inserts
            s = edlib.align(A1, A2[a:e], task='distance')
            distance = s['editDistance']
            # distance = l_distance(A1, A2[a:e])
            iden = round(1 - distance/len(A2[a:e]),3)
            if iden >= idenprev:
                templist.append([a,iden])
                idenprev = iden
            else:
                break
        templist.sort(key=lambda x: x[1]) # sort according to value
        a, iden = templist[-1] # get best value
    else:
        # distance = l_distance(A1, A2)
        s = edlib.align(A1, A2, task='distance')
        distance = s['editDistance']
        iden = round(1 - distance/len(A2),3)
        
    return iden
#==============================================================================
def compl_reverse(self):
    inp  = 'ATCG' # translate table for complement
    outp = 'TAGC'
    complement = ''.maketrans(inp, outp)
    R = (self[::-1]).translate(complement)  # complement reverse
    return R
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
        '_total_outputfig.pdf')
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
    with open(os.path.join(infolder,infile), 'r') as inf: # check the fileformat
        line = inf.readline()
        if line[0] == '>':
            fileformat = 'fasta'
        elif line[0] == '@':
            fileformat = 'fastq'
            
    minlength = args.minlength
    maxlength = args.maxlength
    maxreads = args.maxreads
    ran = args.random
    all_reads = args.all # compare all reads with each other
    global comparelist2, num_seq, comparelist
    comparelist = []
    readlengthlist = []
    ti = 0  # total number of reads in the file
    print('Reading ' + self)
    inputfile = open(os.path.join(infolder,infile), "r")
    for record in SeqIO.parse(inputfile, fileformat):
        ti += 1
        readlengthlist.append(len(record.seq)) # for figure
        # add name, sequence, index number to list
        # index number is needed for resorting to original order later 
        if len(record.seq) >= minlength:
            if maxlength is None:
                comparelist.append([record.id, str(record.seq), '', '']) 
            else:
                if len(record.seq) <= maxlength:
                    comparelist.append([record.id, str(record.seq), '', ''])
    total_num_seq = len(comparelist) # total number of reads passed selection
    
    for i,[w,x,y,z] in enumerate(comparelist):
        comparelist[i] = [w, x, 'u', i] # replace third pos with 'u' (unique) 
                                        # and 4th pos of each list item with 
                                        # index number
        
    comparelist2 = []
    
    if total_num_seq < 1000:
         comparelist2.append(comparelist)
         sentence = '--> Low number of reads, reading all '
    else:
        if ran == True:
            if all_reads == True:
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
            if all_reads == True:
                if maxreads > total_num_seq:
                    print('Comparing all reads (-a, -all) with each other is not '
                      'possible when selecting more reads than available ('
                      + str(maxreads) + ' selected - ' + str(total_num_seq) + 
                      ' available)')
                    sys.exit()
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
        sys.exit()
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
def process_list(self): # make files to do comparisons
    global comparelist2, len_todolist
    nprocesses = args.nprocesses
    # compare 1 with 2,3,4,5,...; compare 2 with 3,4,5,... compare 3 with 4,5,...
    chunk = 500000  # split size of the chunks to feed the multiprocessing queue
    len_todolist = 0
    todoqueue = Queue(maxsize = 2)
    outputfolder = args.outputfolder
    try:  # remove temporary file if exists
        for x in glob.glob(os.path.join(outputfolder, '*.todo')):
            os.remove(x)
        time.sleep(1)
    except FileNotFoundError:
        pass
    
    def queuer(): # function to place items in files for the queue
        global len_todolist
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
                    time.sleep(120)  
                    for dirpath, dirnames, filenames in os.walk(outputfolder):
                        filenames = [i for i in filenames if i.endswith('.todo')]
        todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
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
            process = [Process(target=similarity, args=(todoqueue,)) for x in 
                       range(nprocesses)]
            for p in process:
                p.start()
            for p in process:
                p.join() 
        except KeyboardInterrupt:
            print("Shutting processes down")
           # Optionally try to gracefully shut down the worker processes here.
            p.terminate()
            p.join()
    
    Thread(target = queuer).start()
    Thread(target = feeder).start()
    time.sleep(10)
    c = Thread(target = consumer)
    c.start()
    c.join() # wait until c has finished its work
#==============================================================================
def similarity(todoqueue): # process files for similarity 
    global progress, len_todolist
    try:  # remove temporary file if exists
        os.remove(tempfile)
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
        with open(tempfile, 'a') as f:
            for c in templist:
                f.write(c +'\n')
        templist =[]
        MYLOCK.release()
#==============================================================================
def SSG(tempfile):  #calculate the N6
    # comparable with the N50 in sequence assembly, here used to estimate ssg value
    print('Estimating the ssg value for this dataset')
    t = 0
    totalsimil = 0
    tempdict = {}
    with open(tempfile, 'r') as tf:
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
    global consensusgroup, num_seq 
    templist = []
    tempdict = {}
    t = 0 # items done
    try:
        with open(tempfile, 'r') as tf:
            for line in tf:
                e = line.strip().split(':')[:3]
                if e[1] in tempdict:
                    tempdict[e[1]].append(e)
                    # sort list based on 2nd number (A2) and score)
                    tempdict[e[1]].sort(key=lambda x: (int(x[1]), float(x[2]))) 
                    for i, j in enumerate(tempdict[e[1]][:-1]):
                        # if second index number is the same for the next item
                        if j[1] == tempdict[e[1]][i+1][1]: 
                            # if iden is lower: keep the best
                            if j[2] < tempdict[e[1]][i+1][2]: 
                                j[0] = ''  # mark to remove
                    # remove those lower values from the list 
                    tempdict[e[1]] = [i for i in tempdict[e[1]] if i[0] != ''] 
                else:
                    tempdict[e[1]] = [e]
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
        print('--------------------------------------------------------------------------')
        print('"infile_"compare.tmp and "infile"_comparelist.pickle are missing.')
        print('You have to do process the file without the "--species_only" option first.')
        print('--------------------------------------------------------------------------')
        sys.exit()
    
    grouplist = []  
    # Make groups with sequences with high similarity
    for x in templist:  
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
                                  replace('.fasta', '_') + str(j) + '.group')
        with open(outputfile, 'a') as outputf:
            for y in x:
                outputf.write(str(y) + '\n')
        print('  ' + str(os.path.split(outputfile)[1]) + ' contains ' + 
              str(len(x)) + ' sequences (' + str(round(len(x)*100/num_seq, 2)) + '%)')
    print(str(grouped_seq) + '/' + str(num_seq) + ' sequences assigned in groups (' 
          + str(round(grouped_seq*100/num_seq, 2)) + '%)')

    return grouplist, args 
#==============================================================================
def merge_groups(grouplist):
    a1 = len(grouplist)
    a2 = 0
    print('--> Number of groups before merge: ' + str(a1)) 
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
def comp_consensus_groups(grouplist): # compare consensuses with each other
    global comparelist
    comparelist.sort(key=lambda x: x[3]) # sort list based on index number
                                         # must be done for option '-all' sequences !
    print('-> Merging based on consensus of 50 reads per group') 
    grouplist = [list(i) for i in grouplist]  # make list of the groupsets
    a1 = len(grouplist)
    a2 = 0
    while a1 > a2:
        a1 = len(grouplist)
        position = 0
        templist = []
        for i,x in enumerate(grouplist): # create consensus of 50 seq in each group
            consensuslist = []
            if len(x) > 50: # if number of reads > 50, only take 50 for consensus
                x = random.sample (x, 50)
            # get the sequence that matches the number
            for y in x:
                consensuslist.append(comparelist[int(y)][1]) 
            consensuslist.sort(key=lambda x: len(x)) #sort list based on length seq
            # get all seq in same direction 
            consensuslist2 = consensus_direction(consensuslist) 
            consensus = l_median(consensuslist2)  # create consensuse sequence
            grouplist[i].append(consensus)  # add consensus sequence to the group 
            print(str(round(i/a1*100, 1)) + '% consensuses made', end='\r') 
            # to clear previous line completely   
        print('                                                ', end='\r')  
        y = 0 #position in first range
        z = 0 #position in 2nd range
        for position in range(position, len(grouplist)-1):
            z = y
            for position2 in range(position+1,len(grouplist)):
                z +=1
                A1 = grouplist[position]
                A2 = grouplist[position2]
                if len(A1[-1])*1.08 < len(A2[-1]) or len(A2[-1])*1.08 < len(A1[-1]):
                    pass  #don't compare if length of sequences differ to much
                else:
                    idenlist = []
                    iden = distance(A1[-1],A2[-1])
                    idenlist.append(iden) # add iden to list
                    idenR = distance(A1[-1],compl_reverse(A2[-1]))
                    idenlist.append(idenR) # add idenR to list
                    idenlist.sort(reverse=True) # sort the list
                    iden = idenlist[0] # take the biggest value
                    if iden >= 0.60 : 
                        templist.append([y, z, iden])
            y += 1  
            print(str(round(y/a1*100, 1)) + '% comparisons done', end='\r') 
        for x in templist:  
            # extend both groups with each other,
            grouplist[x[0]].extend(grouplist[x[1]])  
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
    indexes = set()
    outputfolder = args.outputfolder
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
    templist = []
    try:
        with open(tempfile, 'r') as tf:
            for line in tf:
                line = line.strip().split(':')
                if len({line[0], line[1]}.intersection(indexes)) > 0:
                    if float(line[2]) >= similar_species_groups: 
                        templist.append(line)
    except FileNotFoundError:
        pass
    # only keep the best values:
    a = len(templist)
    b = 0
    while a > b:
        a = len(templist)
        #sort list based on 2nd number (A2) and score
        templist.sort(key=lambda x: (int(x[1]), float(x[2])))  
        for i, j in enumerate(templist[:-1]):
            if j[1] == templist[i+1][1]: # if second index number is the same 
                if j[2] < templist[i+1][2]: # if iden is lower
                    j[0] = ''                # mark to remove
        templist = [j for j in templist if j[0] != ''] # remove from list  
        b = len(templist)
    
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
    grouplist = [list(set(i)) for i in grouplist if len(i) > 2] 
    a2 = len(grouplist)
    print('--> Number of groups after removing groups with less than 3 sequences: '
          + str(a2)) 
    
    print('----> Making consensus for each group')
    
    comparelist.sort(key=lambda x: x[3]) # sort list based on index number
                                         # must be done for 'all' sequences !
    for i,x in enumerate(grouplist):
        consensuslist = []
        if len(x) > 500: # if number of reads > 500, only take 500 for consensus
            x = random.sample (x, 500)
        for y in x:
            # get the sequence that matches the number
            consensuslist.append(comparelist[int(y)][1]) 
        # get all seq in same direction
        consensuslist2 = consensus_direction(consensuslist) 
        consensus = l_median(consensuslist2)  # create consensuse sequence
        grouplist[i].append(consensus)  # add consensus sequence to the group

    return indexes, grouplist
#==============================================================================
def filter_seq(group_filename, grouplist, indexes):
    # filter sequences: put the sequences with high similarity in separate files.
    # Sequences of the same species with the same gene should be in one file.
    with open(os.path.join(infolder, infile), 'r') as inf: # check the fileformat
        line = inf.readline()
        if line[0] == '>':
            fileformat = 'fasta'
        elif line[0] == '@':
            fileformat = 'fastq'
    outputfolder = args.outputfolder    
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
        record_dict = SeqIO.index(os.path.join(infolder, infile), 'fastq')  
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
            with open(outputfile, 'a') as outputf:
                x = str(seq)
                outputf.write('>' + str(index)  + '\n' + x + '\n')
                if fq == True:
                    with open(outputfilefq, 'a') as writer: 
                        SeqIO.write(record_dict[rec_id], writer, 'fastq')
    MYLOCK.release()
    grouped_seq = 0 # number of sequences in grouplist  
    for x in grouplist:
        for y in x:
            if y.isdigit(): # list contains consensus sequence, don't count those
                grouped_seq += 1
                
    consensusfilename = os.path.join(outputfolder, group_filename).replace(
        '.group', '_consensussequences.fasta') # group consensusfile
    # total consensusfile
    consensusfile = os.path.join(outputfolder, infile).replace('.fasta', 
    '_consensussequences.fasta').replace('.fastq', '_consensussequences.fasta')
    
    try:  # remove  file if exists
        os.remove(os.path.join(outputfolder, consensusfilename))
    except FileNotFoundError:
        pass 

    MYLOCK.acquire()       
    for j, x in enumerate(grouplist):
        outputfile = os.path.join(outputfolder, group_filename).replace(
            '.group', '_') + str(j) + '.fasta'
        consensusname = group_filename.replace('.group', '_') + str(j) 
        
        t = 0
        for y in x:
            if y.isalpha(): # if it is a sequence
                with open(outputfile, 'a') as outputf:
                    outputf.write('>consensus' + '\n' + y + '\n')
                with open(consensusfilename, 'a') as outputf:
                    outputf.write('>consensus_' + str(consensusname) + '(' +
                                  str(len(x)-1) + ')\n' + y + '\n')
                with open(consensusfile, 'a') as outputf:
                    outputf.write('>consensus_' + str(consensusname) + '(' + 
                                  str(len(x)-1) + ')\n' + y + '\n')

            elif y.isdigit():
                t += 1 # count number of sequences in group
        with open(os.path.join(outputfolder, 'results.txt'), 'a') as rf:
            try:
                rf.write('--> ' + str(os.path.split(outputfile)[1]) + ' contains ' 
                         + str(t) + ' sequences (' + 
                         str(round((len(x)-1)*100/len(indexes), 2)) + '%)\n')
            except ZeroDivisionError:
                pass
    with open(os.path.join(outputfolder,'results.txt'), 'a') as rf:
        try:
            rf.write(str(grouped_seq) + '/' + str(len(indexes)) + 
                     ' sequences assigned in groups for ' + group_filename + 
                     ' (' + str(round(grouped_seq*100/len(indexes), 2)) + '%)\n')
        except ZeroDivisionError:
            pass
    MYLOCK.release()
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
                if len(todolist) == 500000: # save in chuncks to save memory
                    todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
                    with open(todofilename, 'wb') as wf:
                        pickle.dump(todolist, wf)
                        todolist = []
                        k += 1
    todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
    with open(todofilename, 'wb') as wf:
        pickle.dump(todolist, wf)

    print(group_filename + '----> ' + str(l) + ' comparisons to calculate')

    try:  # remove temporary file if exists
        os.remove(os.path.join(outputfolder, group_tempfile))
    except FileNotFoundError:
        pass    

    for dirpath, dirnames, filenames in os.walk(outputfolder):
        filenames = [i for i in filenames if i.endswith('.todo')]
        filenames.sort(key=lambda x: os.path.getmtime(os.path.join(outputfolder, x)))
        for name in filenames:
            print('...processing: ' + name)
            with open(os.path.join(outputfolder, name), 'rb') as rf:
                todolist = pickle.load(rf)
                if len(todolist) >= nprocesses:
                    chunk = len(todolist)//nprocesses
                else:
                    chunk = len(todolist)
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
                    process = [Process(target=similarity_species, args=(todoqueue, 
                            outputfolder, group_tempfile,)) for x in range(nprocesses)]
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
def similarity_species(todoqueue, outputfolder, group_tempfile):
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
    templist = []
    templist2 = []
    min_similar = args.similar_species/100
    outputfolder = args.outputfolder
    group_tempfile = os.path.join(outputfolder, group_filename).replace('.group', '.tmp')
    try:
        with open(group_tempfile, 'r') as tf:
            temp = tf.readlines()
        for line in temp:
            templist.append(line.strip().split(':'))
    except FileNotFoundError:
        pass
    # only keep the best values:
    a = len(templist)
    b = 0
    while a > b:
        a = len(templist)
        #sort list based on 1st number and score
        templist.sort(key=lambda x: (int(x[0]), float(x[2])))  
        for i, j in enumerate(templist[:-1]):
            if j[0] == templist[i+1][0]: # if second index number is the same 
                if j[2] <= templist[i+1][2]: # if iden is lower or equal
                    j[0] = ''                # mark to remove
        templist = [j for j in templist if j[0] != ''] # remove from list  
        b = len(templist)
        
    templist2 = [i for i in templist if float(i[2]) >= similar] #only keep those 
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
   
        for i,x in enumerate(grouplist):
            if x[-1].isdigit(): # if the last item is a number, sequence has been added
                consensuslist = []
                if len(x) > 100: # if number of reads > 100, only take 100 for consensus
                    x = random.sample(x, 100)
                for y in x:
                    # get the sequence that matches the number
                    consensuslist.append(comparelist[int(y)][1]) 
                # get all seq in same direction
                consensuslist2 = consensus_direction(consensuslist) 
                consensus = l_median(consensuslist2)  # create consensus sequence
                grouplist[i].append(consensus)  # add consensus sequence to the group
        
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
def compare_consensus(grouplist):
    # compare consensusses with each other
    similar_consensus = args.similar_consensus/100    
    a1 = len(grouplist)
    a2 = 0
    while a1 > a2:
        a1 = len(grouplist)
        position = 0
        grouplist2 = grouplist.copy() # need copy to delete and extend items
        y = 0 #position in first range
        z = 0 #position in 2nd range
        for position in range(position, len(grouplist)-1):
            z = y
            for position2 in range(position+1,len(grouplist)):
                z +=1
                A1 = grouplist[position]
                A2 = grouplist[position2]
                if len(A1[-1])*1.08 < len(A2[-1]) or len(A2[-1])*1.08 < len(A1[-1]):
                    pass  #don't compare if length of sequences differ to much
                else:
                    idenlist = []
                    iden = distance(A1[-1],A2[-1])
                    idenlist.append(iden) # add iden to list
                    idenR = distance(A1[-1],compl_reverse(A2[-1]))
                    idenlist.append(idenR) # add idenR to list
                    idenlist.sort(reverse=True) # sort the list
                    iden = idenlist[0] # take the biggest value
                    if iden >= similar_consensus : 
                        grouplist2[position].extend(A2)
            y += 1  
                     
        grouplist = merge_groups(grouplist2)
    
        for i,x in enumerate(grouplist):  # remove consensus from groups
            grouplist[i] = [y for y in x if y.isdigit()]
        
        for i,x in enumerate(grouplist):
            consensuslist = []
            if len(x) > 200: # if number of reads > 200, only take 200 for consensus
                x = random.sample(x, 200)
            for y in x:
                # get the sequence that matches the number
                consensuslist.append(comparelist[int(y)][1]) 
            # get all seq in same direction
            consensuslist2 = consensus_direction(consensuslist) 
            consensus = l_median(consensuslist2)  # create consensuse sequence
            grouplist[i].append(consensus)  # add consensus sequence to the group
        a2 = len(grouplist)
        
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
    process_consensuslist(indexes, grouplist, group_filename)                
    comparelist2, grouplist, templist, templist2 = update_groups(group_filename, 
                                                                 grouplist)
    k += 1
    while similar > min_similar: 
        while k <= 2:
            if len(templist2) > 0:
                process_consensuslist(indexes, grouplist, group_filename)
                comparelist2, grouplist, templist, templist2 = update_groups(
                    group_filename, grouplist)
                if len(templist2) > 0 and k < 2:
                    # compare consensusses with each other
                    grouplist = compare_consensus(grouplist)  
                k += 1
            else: 
                k = 3
        else:
            k = 0
            similar = round(similar - 0.01, 2) 
            print(group_filename + '----> similarity = ' + str(similar))
            if len(templist) == 0:
                process_consensuslist(indexes, grouplist, group_filename)
                comparelist2, grouplist, templist, templist2 = update_groups(
                    group_filename, grouplist)
                # grouplist = compare_consensus(grouplist) # compare consensusses with each other 
                k += 1
            else: 
                comparelist2, grouplist, templist, templist2 = update_groups(
                    group_filename, grouplist)
                k += 1
    grouplist = compare_consensus(grouplist) # compare consensusses with each other 

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
        filename = infile.replace('.fastq', '_*').replace('.fasta', '_*')
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
        process_list(comparelist2) 
        
        comparelist2 = list(set(([tuple(x) for sublist in comparelist2 for x in 
                                  sublist]))) # make list out of list with sublists
        comparelist2 = [list(x) for x in comparelist2]
        comparelist2.sort(key=lambda x: x[3])
        
        with open(saved_comparelist, 'wb') as wf:
            pickle.dump(comparelist, wf)
            pickle.dump(comparelist2, wf)
        with open('readme.txt', 'w') as rm:
            print('This folder contains 2 temporary files:', file = rm)
            print('- "infile"_compare.tmp = text file with similarity between '
                  'numbered sequences.', file = rm)
            print('- "infile"_comparelist.pickle = binary file with all sequence'
                  ' info in python', file = rm)
            print('  pickle format.  ', file = rm)
            print('', file = rm)
            print('If you want to play with the "--similar_species_groups" or '
                  '"--similar_species"', file = rm)
            print('parameters, you can save time if these files are kept here.', 
                  file = rm)
            print('The "-so --species_only" option is used to do this.', file = rm)
            print('', file = rm)
            print('Otherwise, these 2 files can be removed.', file = rm)
#==============================================================================
def sort_groups(): # read the gene groups and sort sequences to species level
    global comparelist, comparelist2, resultlist, num_seq, saved_comparelist
    outputfolder = args.outputfolder
    try:
        with open(saved_comparelist, 'rb') as rf:
            comparelist = pickle.load(rf)
            comparelist2 = pickle.load(rf)
        num_seq = len(comparelist2)
    except FileNotFoundError:
        pass 
    
    consensusfile = infile.replace('.fasta', '_consensussequences.fasta').replace(
        '.fastq', '_consensussequences.fasta') # total consensusfile
    try:  # remove  file if exists
        os.remove(os.path.join(outputfolder, consensusfile))
    except FileNotFoundError:
        pass   
    
    update_list(tempfile)
               
    for dirpath, dirnames, filenames in os.walk(outputfolder):
        for name in filenames:
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
        check_version(version)
        if args.similar_species_groups is None:
            args.similar_species_groups = 'Estimate'
        save_arguments() # write all settings in the results.txt file
        infolder_file = args.input
        infolder, infile = os.path.split(os.path.realpath(infolder_file))
        tempfile = infile.replace('.fastq', '_compare.tmp').replace('.fasta', 
                                                                    '_compare.tmp')
        saved_comparelist = infile.replace('.fastq', '_comparelist.pickle').replace(
            '.fasta', '_comparelist.pickle')
        if args.species_only == True:
            sort_groups()
        else:
            sort_genes()
            sort_groups()
    except KeyboardInterrupt:
        sys.exit()