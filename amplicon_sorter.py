# -*- coding: utf-8 -*-
"""

Python 3.5.1
(March 2019)

File version 2020.03.11 

@author: Andy Vierstraete

compare sequences from a MinION run and sort them based on similarity

For c implementation of Levenshtein:
    https://pypi.org/project/python-Levenshtein/
    sudo apt-get install python3-dev
    (needed: python3-setuptools, python3-pip, python3-wheel)
    python3 -m pip install python-Levenshtein

"""
from Bio import SeqIO
import multiprocessing
from multiprocessing import Process, Lock, Queue
from threading import Thread
from Levenshtein import distance as l_distance
from Levenshtein import median as l_median
import os
import random
import time
import glob
import pickle
import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

global tempfile, infile, num_seq
global chunk
tempfile = 'compare1.tmp'
#tempfile2 = 'compare2.tmp'
chunk = 500000   # split size of the chunks to feed the multiprocessing queue
#similar = 0.5
#infile = 'rolish7_trim.fastq'




#==============================================================================
#def levenshtein_distance(s1, s2):  # function to compare 2 sequences
#    """
#    Python version of Levenshtein distance for compatability.  Taken from:
#    http://code.activestate.com/recipes/576874-levenshtein-distance/
#    Little adjustments for Python 3.
#    """
#    l1 = len(s1)
#    l2 = len(s2)
#
#    matrix = [range(l1 + 1)] * (l2 + 1)
#    for zz in range(l2 + 1):
#        matrix[zz] = list(range(zz, zz + l1 + 1))
#    for zz in range(0, l2):
#        for sz in list(range(0, l1)):
#            if s1[sz] == s2[zz]:
#                matrix[zz + 1][sz + 1] = min(matrix[zz + 1][sz] + 1,
#                                    matrix[zz][sz + 1] + 1, matrix[zz][sz])
#            else:
#                matrix[zz + 1][sz + 1] = min(matrix[zz + 1][sz] + 1,
#                                    matrix[zz][sz + 1] + 1, matrix[zz][sz] + 1)
#    return matrix[l2][l1]
#==============================================================================
def get_arguments():
    
    def range_limited_float_type(arg):
        """ Type function for argparse - a float within some predefined bounds """
        try:
            f = float(arg)
        except ValueError:    
            raise argparse.ArgumentTypeError("Must be a floating point number")
        if f < 50 or f > 100:
            raise argparse.ArgumentTypeError("Argument must be > " + str(50.0) + " and < " + str(100.0))
        return f

    def valid_file(param):  # check if input file ends on .fastq or .fasta
        base, ext = os.path.splitext(param)
        if ext.lower() not in ('.fasta', '.fastq'): 
            raise argparse.ArgumentTypeError('File extension must be .fastq or .fasta') 
        return param
        
    parser = argparse.ArgumentParser(description='AmpliconSorter: Sort amplicons based on identity ' 
                                     'and saves them in different files including the consensus.' )
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
    parser.add_argument('-sg', '--similar_genes', type = range_limited_float_type, required=False, default=50.0 ,
                        help='Similarity to sort genes in groups (value between 50 and 100). Default=50.0')
    parser.add_argument('-ss', '--similar_species', type = range_limited_float_type, required=False, default=85.0 ,
                        help='Similarity to sort to species level (value between 50 and 100). Default=85.0')
    parser.add_argument('-sfq', '--save_fastq', action = 'store_true',
                        help='Save the results also in fastq files (fastq files will not contain the consensus sequence)')
    parser.add_argument('-ra', '--random', action = 'store_true',
                        help='Takes random reads from the inputfile if --maxreads is lower than total number of reads that passed criteria')

    args = parser.parse_args()
    return args
#==============================================================================    
def distance(A1,A2):  # calculate the similarity of 2 sequences
    if len(A1)*1.05 < len(A2): # if A1 is much shorter than A2, it influences the similarity
        idenprev = 0
        templist = []
        for a in range(0, len(A2) - len(A1),2): 
            e = int(len(A1)*1.015) # leave room for inserts
            distance = l_distance(A1, A2[a:e])
            iden = round(1 - distance/len(A2[a:e]),3)
            if iden >= idenprev:
                templist.append([a,iden])
                idenprev = iden
            else:
                break
        templist.sort(key=lambda x: x[1]) # sort taccording to value
        a, iden = templist[-1] # get besst value
    else:
        distance = l_distance(A1, A2)
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
def figure(readlengthlist):
    bases = sum(readlengthlist)
    length_N50 = N50(readlengthlist, bases) 
    try:
        max_readlength = sorted(readlengthlist)[-1]
    except IndexError:
        pass
    plt.figure(1, figsize=[5,5])
    ax = plt.subplot(2,1,1)
    figname = infile.replace('.fastq', '_total_outputfig.pdf').replace('.fasta', '_total_outputfig.pdf')
    plt.ylabel('Number of reads')
    plt.title('Read length histogram') 
    plt.text(0.95, 0.55, 'Total yield (Gb): ' + str(round((bases/1000000000),2)) + 
             '\nNumber of reads: ' +
             '{:,}'.format(len(readlengthlist)) + 
             '\nMax readlength: ' + str(round(max_readlength/1000,1)) + ' Kb' \
             '\nN50 = ' + str(round((length_N50/1000),1)) + ' Kb',   
             horizontalalignment='right', transform=plt.gca().transAxes)
    
    plt.hist(readlengthlist, bins='auto', color='green')
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:,.1F}'.format(x/1000))) # divide read length by 1000
    plt.savefig(figname, format='pdf', dpi=300)

    ax2 = plt.subplot(2,1,2)
    plt.ylabel('Log Number of reads')
    plt.xlabel('Read length (Kb)')
   
    plt.hist(readlengthlist, log=True, bins='auto', color='green')
    ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:,.1F}'.format(x/1000))) # divide read length by 1000
#    plt.legend(loc='center right')
    plt.tight_layout()
    plt.savefig(figname, format='pdf', dpi=300)
    plt.clf() # clear the figure            

    print('Saved "' + figname + '" as a Read Length Histogram.')

#==============================================================================
def read_file(self):
    ''' 
    options to include:
    -read all files from a folder
    -input can be fastq or fasta, autodetect file type
    '''
    global fileformat
    with open(infile, 'r') as inf: # check the fileformat
        line = inf.readline()
        if line[0] == '>':
            fileformat = 'fasta'
        elif line[0] == '@':
            fileformat = 'fastq'
            
    minlength = args.minlength
    maxlength = args.maxlength
    maxreads = args.maxreads
    ran = args.random
    global comparelist, num_seq
    comparelist = []
    readlengthlist = []
    ti = 0  # total number of reads in the file
    print('Reading ' + self)
    with open (logfile, 'w') as lf:
        print('Reading ' + self, file = lf)
    inputfile = open(self, "r")
    for record in SeqIO.parse(inputfile, fileformat):
        ti += 1
        readlengthlist.append(len(record.seq))
        # add name, sequence, index number to list
        # index number is needed for resorting to original order later (filter_seq) to get
        # higher speed when searching for item in list
        if len(record.seq) >= minlength:
            if maxlength is None:
                comparelist.append([record.id, str(record.seq), '', '']) #new
            else:
                if len(record.seq) <= maxlength:
                    comparelist.append([record.id, str(record.seq), '', ''])
    total_num_seq = len(comparelist) # total number of reads passed selection
    if ran == True:
        comparelist = random.sample(comparelist, maxreads)
        sentence = '--> Reading random '
    else:
        comparelist = comparelist[0:maxreads] 
        sentence = '--> Reading '
    inputfile.close()
    print(self + ' contains ' + str(ti) + ' reads.')
    figure(readlengthlist) # make a read length histogram from the inputfile
    comparelist.sort(key=lambda x: len(x[1])) #sort list based on length seq
    for i,[w,x,y,z] in enumerate(comparelist):
        comparelist[i] = [w, x, 'u', i] # replace third pos with 'u' (unique) and 4th pos of each 
                                        # list item with index number
    num_seq = len(comparelist) # number of sample reads
    if maxlength is None:
        print(sentence + str(num_seq) + ' out of ' + str(total_num_seq) + 
              ' sequences longer than ' + str(minlength) + 'bp')
    else:
        print(sentence + str(num_seq) + ' out of ' + str(total_num_seq) + 
              ' sequences between ' + str(minlength) + ' and ' + str(maxlength) + 'bp')
#    with open (logfile, 'a') as lf:
#        print('--> File contains ' + str(num_seq) + ' sequences' , file = lf)
    comparelist = comparelist[:]  # for debugging with smaller list
    return comparelist, num_seq
#==============================================================================
def process_list(self):
    global comparelist, len_todolist, progressfile
    nprocesses = args.nprocesses
    # compare 1 with 2,3,4,5,...; compare 2 with 3,4,5,... compare 3 with 4,5,...

    len_todolist = 0
    todoqueue = Queue(maxsize = 2)
    progressfile = 'progress.tmp'
    try:  # remove temporary file if exists
        os.remove(progressfile)
    except FileNotFoundError:
        pass
    try:  # remove temporary file if exists
        for x in glob.glob('*.todo'):
            os.remove(x)
        time.sleep(1)
    except FileNotFoundError:
        pass
    
    def queuer(): # function to place items in files for the queue
        global len_todolist
        
#        print('start queueing')
        position = 0
        y = 0 #position in first range
        z = 0 #position in 2nd range
        l = 0 # length of todolist
        k = 0 # number of todofiles
        tl = 0 #total length todolist
        todolist = []
        for position in range(position, len(self)-1):
            z = y
            for position2 in range(position+1,len(self)):
                z +=1
                A1 = self[position]
                A2 = self[position2]
                if len(A1[1])*1.3 < len(A2[1]):
                    pass  #don't compare if length of sequences differ to much
                else:
                    todolist.append([A1,A2])
                    l += 1
                    tl += 1
                    if l == chunk:
                        todofilename = 'file_' + str(k) + '.todo'
                        with open(todofilename, 'wb') as wf:
                            pickle.dump(todolist, wf)
                        todolist = []
                        k += 1
                        l = 0
            y += 1
            todofilename = 'file_' + str(k) + '.todo'
            with open(todofilename, 'wb', buffering=0) as wf:
                pickle.dump(todolist, wf)
            for dirpath, dirnames, filenames in os.walk(os.getcwd()):
                filenames = [i for i in filenames if i.endswith('.todo')]
            while len(filenames) > 20: # limit numbers of temporary files
                time.sleep(120)  
                for dirpath, dirnames, filenames in os.walk(os.getcwd()):
                    filenames = [i for i in filenames if i.endswith('.todo')]
        len_todolist = tl
       
        return len_todolist
     
    def feeder(): # function to feed the queue in parts to save memory
        global len_todolist
        while len_todolist == 0:  #start feeding when it is still making todolist
#            print('start feeding')
            time.sleep(30)
            for dirpath, dirnames, filenames in os.walk(os.getcwd()):
                filenames = [i for i in filenames if i.endswith('.todo')]
                filenames.sort(key = os.path.getmtime) #sort on modification time
                filenames = filenames[:-1] # don't include last file, possible still being written to disk
                for name in filenames[:3]: # only take a few files in memory to get out of while loop asap
                    print('processing: ' + name)
                    with open(name, 'rb') as rf:
                        sublist = pickle.load(rf)
                        todoqueue.put(sublist, block=True) 
                        time.sleep(2)
                        os.remove(name)
        else: # feed all the rest when finished making todolist
            for dirpath, dirnames, filenames in os.walk(os.getcwd()):
                filenames = [i for i in filenames if i.endswith('.todo')]
                filenames.sort(key = os.path.getmtime)
                for name in filenames:
                    print('processing: ' + name)
                    with open(name, 'rb') as rf:
                        sublist = pickle.load(rf)
                        todoqueue.put(sublist, block=True)
                        time.sleep(2)
                        os.remove(name)
            for i in range(nprocesses): # put 'STOP' at the end of the queue for every process
                todoqueue.put("STOP")    
                        
    def consumer(): # function to consume the queue  
#        print('start consuming')
        try:
            process = [Process(target=similarity, args=(todoqueue,)) for x in range(nprocesses)]
#            for p in process:
#                #ask the processes to stop when all files are handled
#                #"STOP" is at the very end of queue
#                todoqueue.put("STOP")Thread(target = 
            for p in process:
                p.start()
            for p in process:
                p.join() 
            #for i in range(nprocesses): # put 'STOP' ant the end of the queue for every process
       
        except KeyboardInterrupt:
            print("Shutting processes down")
           # Optionally try to gracefully shut down the worker processes here.
            p.terminate()
            p.join()
    
   
    Thread(target = queuer).start()
    Thread(target = feeder).start()
    time.sleep(60)
    c = Thread(target = consumer)
    c.start()
    c.join() # wait until c has finished its work
    
#==============================================================================
def similarity(todoqueue):
    global progress, chunk, len_todolist
    try:  # remove temporary file if exists
        os.remove(tempfile)
    except FileNotFoundError:
        pass
    similarg = args.similar_genes/100
    templist = []
    q = 0 # number of comparisons done
    MYLOCK = Lock()
    for X in iter(todoqueue.get, 'STOP'):    #do stuff until infile.get returns "STOP"
        for A1, A2 in X:
            q += 1
            try:
                b = multiprocessing.current_process()
                # list items are: rec_id, seq, index
                iden = distance(A1[1],A2[1])
                if iden >= similarg: # limit of being the same gene for nanopore data ?
#                    print('Similarity is ' + str(iden))
                    templist.append((str(A1[3]) + ':' + str(A2[3]) + ':' + str(iden)))
                elif iden < similarg:
                    iden = distance(A1[1],compl_reverse(A2[1]))
                    if iden >= similarg:
#                        print(str(A1[3]) + ':' + str(A2[3]) + ':' + str(iden) + ':' + ' reverse')
                        templist.append((str(A1[3]) + ':' + str(A2[3]) + ':' + str(iden) + ':' + 'reverse'))
            except KeyboardInterrupt:
                print("Shutting process 1 down")
                b.terminate()
        MYLOCK.acquire()  # save the list every chunk similarities
        with open(tempfile, 'a') as f:
            for c in templist:
                f.write(c +'\n')
        templist =[]

        q = 0
        MYLOCK.release()

#==============================================================================
def update_list(tempfile):
    similarg = args.similar_genes/100
    print('Filtering compared sequences for best hits and create groups')
    with open (logfile, 'a') as lf:
        print('Filtering compared sequences for best hits and create groups', file = lf)
    global comparelist, consensusgroup, num_seq #, grouplist
    templist = []
    t = 0 # items in templist
    try:
        with open(tempfile, 'r') as tf:
            for line in tf:
                templist.append(line.strip().split(':'))
                t += 1
                if t == 1000000: # clean up the list from time to time to save memory
                    t = 0
                    templist.sort(key=lambda x: (str(x[1]), str(x[2]))) #sort list based on 2nd number (A2) and score 
                    for i, j in enumerate(templist[:-1]):
                        if j[1] == templist[i+1][1]: # if second index number is the same
                            if j[2] <= templist[i+1][2]: # if iden is lower or equal
                #                if float(j[2]) < similar:
                                j[2] = ''                # mark to remove
                    templist = [i for i in templist if i[2] != ''] #only keep those with highest iden  
        templist.sort(key=lambda x: (str(x[1]), str(x[2]))) #sort list based on 2nd number (A2) and score 
        for i, j in enumerate(templist[:-1]):
            if j[1] == templist[i+1][1]: # if second index number is the same
                if j[2] <= templist[i+1][2]: # if iden is lower or equal
    #                if float(j[2]) < similar:
                    j[2] = ''                # mark to remove
        templist = [i for i in templist if i[2] != ''] #only keep those with highest iden  
        templist.sort(key=lambda x: (float(x[2]),int(x[0])),reverse=True) #sort list based on score and index number
        with open (logfile, 'a') as lf:
            print('\n*************************************************************', file = lf)
            print('List with similarities > ' + str(similarg) , file = lf)
            print('*************************************************************',file = lf)
            for x in templist:
                print(x, file = lf)
    except FileNotFoundError:
        pass
    grouplist = []  
    #..........................................................................
    # Make groups with sequences with high similarity
    #..........................................................................
    for j, x in enumerate(templist):   
        for s in grouplist: 
#            if x[1] in s :
#                s.append(x[0]) # extend([x[0], x[1]])
#                templist[j] = [] # remove when added
#                break    
            if x[1] in s or x[0] in s:
                s.extend([x[0], x[1]])
                templist[j] = [] # remove when added
                break    
        else:
            grouplist.append([x[0], x[1]])
 
#    a1 = len(grouplist)
#    a2 = 0
#    print('--> Number of groups before merge: ' + str(a1)) 
#    while a1 > a2:
#        a1 = len(grouplist)
#        y = 0 #position in first range
#        z = 0 #position in 2nd range
#        position = 0
#        for position in range(position, len(grouplist)-1):
#            z = y
#            for position2 in range(position+1,len(grouplist)):
#                z +=1
#                A1 = grouplist[position]
#                A2 = grouplist[position2]
#    #            print(str(y) + ' - ' + str(z) + ' : ' + str((set(A1)).intersection((set(A2)))))
#                if len((set(A1)).intersection((set(A2)))) > 0: # check if numbers occur in other groups
#                    grouplist[position].extend(A2)
#                    grouplist[position2] = [] # mark for removal
#            y += 1 
#        grouplist = [list(set(i)) for i in grouplist if len(i) > 0]  # remove empty sublists
#        a2 = len(grouplist)
#        print('--> Number of groups after merge: ' + str(a2)) 
        
    grouplist = merge_groups(grouplist)
    print('comparing consensuses in the groups')
    grouplist = comp_consensus_groups(grouplist)
        
    grouped_seq = 0 # number of sequences in grouplist  
    for x in grouplist:
        grouped_seq += len(x)  
        
    for j, x in enumerate(grouplist):
        outputfile = infile.replace('.fastq', '_').replace('.fasta', '_') + str(j) + '.group'
        with open(outputfile, 'a') as outputf:
            for y in x:
                outputf.write(str(y) + '\n')
        print('  ' + outputfile + ' contains ' + str(len(x)) + ' sequences (' + str(round(len(x)*100/num_seq, 2)) + '%)')
    print(str(grouped_seq) + '/' + str(num_seq) + ' sequences assigned in groups (' + str(round(grouped_seq*100/num_seq, 2)) + '%)')

    # check which sequences do not belong to a group    
    comparelist.sort(key=lambda x: x[3]) #sort list based on index number
    for o, [rec_id, seq, scores, index] in enumerate(comparelist):
        for i, n in enumerate(grouplist):
            if str(index) in n:  # if it belongs to a group
                comparelist[o][2] = i
    comparelist2 = [i for i in comparelist if i[2] == 'u'] # only keep those that not belong to group
    if len(comparelist2) > 0:
        outputfile = infile.replace('.fastq', '_unique.group').replace('.fasta', '_unique.group')
        with open(outputfile, 'a') as outputf:
            for x in comparelist2:
                outputf.write(str(x[3]) + '\n')
        print('  ' + outputfile + ' contains ' + str(len(comparelist2)) + ' sequences (' + str(round(len(comparelist2)*100/num_seq, 2)) + '%)')
    
    # reset the comparelist
    for i,[w,x,y,z] in enumerate(comparelist):
        comparelist[i] = [w, x, 'u', z] # replace third pos with 'u' (unique) 
#    return grouplist, comparelist#, consensusgroup  
#==============================================================================
def merge_groups(grouplist):
    a1 = len(grouplist)
    a2 = 0
    print('--> Number of groups before merge: ' + str(a1)) 
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
    #            print(str(y) + ' - ' + str(z) + ' : ' + str((set(A1)).intersection((set(A2)))))
                if len((set(A1)).intersection((set(A2)))) > 0: # check if numbers occur in other groups
                    grouplist[position].extend(A2)
                    grouplist[position2] = [] # mark for removal
            y += 1 
        grouplist = [list(set(i)) for i in grouplist if len(i) > 0]  # remove empty sublists
        a2 = len(grouplist)
    print('--> Number of groups after merge: ' + str(a2)) 
    return grouplist
#==============================================================================
def comp_consensus_groups(grouplist):
    position = 0
    grouplist2 = grouplist.copy() # need copy to delete and extend items
    
    for i,x in enumerate(grouplist): # create consensus of 20 seq in each group
        consensuslist = []
        for y in x[0:20]:
            consensuslist.append(comparelist[int(y)][1]) # get the sequence that matches the number
        consensuslist2 = consensus_direction(consensuslist) # get all seq in same direction
        consensus = l_median(consensuslist2)  # create consensuse sequence
        grouplist[i].append(consensus)  # add consensus sequence to the group  
        
    y = 0 #position in first range
    z = 0 #position in 2nd range
    for position in range(position, len(grouplist)-1):
        z = y
        for position2 in range(position+1,len(grouplist)):
            z +=1
            A1 = grouplist[position]
            A2 = grouplist[position2]
            if len(A1[-1])*1.3 < len(A2[-1]) or len(A2[-1])*1.3 < len(A1[-1]):
                pass  #don't compare if length of sequences differ to much
            else:
                idenlist = []
                iden = distance(A1[-1],A2[-1])
                idenlist.append(iden) # add iden to list
                idenR = distance(A1[-1],compl_reverse(A2[-1]))
                idenlist.append(idenR) # add idenR to list
                idenlist.sort(reverse=True) # sort the list
                iden = idenlist[0] # take the biggest value
                if iden >= 0.50 : 
                    grouplist2[position].extend(A2)
#                    with open('consensusseq.txt', 'a') as cs:
#                        print('>1\n' + str(A1[-1]) + '\n>2\n' + str(A2[-1]) + '\n' + 
#                          str(iden) +' '+ str(len(A1[-1])) +' '+ str(len(A2[-1])), file = cs)

        y += 1  
                 
    grouplist = merge_groups(grouplist2)

    for i,x in enumerate(grouplist):  # remove consensus from groups
        grouplist[i] = [y for y in x if y.isdigit()]
    
    return grouplist
#==============================================================================
def read_indexes(group_filename): # read index numbers from the the input file
    global similar, indexes, comparelist #, grouplist
    indexes = set()
    try:
        with open(group_filename, 'r') as gf:
            temp = gf.readlines()
            for line in temp:
                indexes.add(line.strip())
        print('reading ' + group_filename + ' containing ' + str(len(indexes)) + ' sequences.')
    except FileNotFoundError:
        pass
    templist = []
    try:
        with open(tempfile, 'r') as tf:
            for line in tf:
                line = line.strip().split(':')
                if line[0] in indexes and line[1] in indexes:
                    if float(line[2]) >= 0.90:
                        templist.append(line)
    except FileNotFoundError:
        pass
    
    grouplist = []  
    templist.sort(key=lambda x: (float(x[2]),int(x[0])),reverse=True) #sort list based on score and index number
    #..........................................................................
    # Make groups with sequences with high similarity
    #..........................................................................
    while similar >= 0.91:
#        print('similar = ' + str(round(similar, 2)) + ' - # sequences: ' + str(len(templist)))
        for j, x in enumerate(templist):  
            if float(x[2]) >= similar:
                for s in grouplist: 
                    if x[1] in s:
                        s.append(x[0]) # extend([x[0], x[1]])
                        templist[j] = '' # remove when added
                        break    
                else:
                    grouplist.append([x[0], x[1]])
                    templist[j] = '' # remove when added
        templist = [i for i in templist if i != ''] # clean up templist from added sequences
        similar = similar - 0.01
    
    grouplist = merge_groups(grouplist)
    
    print('----> Making consensus for each group')

    for i,x in enumerate(grouplist):
        consensuslist = []
        for y in x:
            consensuslist.append(comparelist[int(y)][1]) # get the sequence that matches the number
        consensuslist2 = consensus_direction(consensuslist) # get all seq in same direction
        consensus = l_median(consensuslist2)  # create consensuse sequence
        grouplist[i].append(consensus)  # add consensus sequence to the group
    
        
#    comparelist.sort(key=lambda x: x[3]) #sort list based on index number
#    for o, [rec_id, seq, scores, index] in enumerate(comparelist):
#        for i, n in enumerate(grouplist):
#            if str(index) in n:  # it belongs to a group
#                comparelist[o][2] = i
    
 
    return indexes, grouplist #, consensusgroup  ,comparelist

#==============================================================================
def filter_seq(group_filename, grouplist, indexes):
    # filter sequences: put the sequences with high similarity in separate files.
    # Sequences of the same species with the same gene should be in one file.
    global fileformat
    if fileformat == 'fasta': # if the inputfile was fasta, it is not possible to 
        fq = False            # save results in fastq format
    else:
        fq = args.save_fastq # check if it needs to be saved in fastq format
    MYLOCK = Lock()
    print('Writing sequences with high similarity in separate files')
    with open (logfile, 'a') as lf:
        print('\nWriting sequences with high similarity in separate files', file = lf)
    global comparelist #, indexes #, grouplist
    MYLOCK.acquire()
    if fq == True:
        record_dict = SeqIO.index(infile, 'fastq')  # index the input fastq file
    for rec_id, seq, scores, index in comparelist:
        if str(index) in indexes:
            if scores == 'u':  # sequences that have no similarity with others
                 outputfile = group_filename.replace('.group', '_') + 'unique.fasta' # unique sequences
                 outputfilefq = group_filename.replace('.group', '_') + 'unique.fastq' # unique sequences
            else:
                 outputfile = group_filename.replace('.group', '_') + str(scores) + '.fasta'
                 outputfilefq = group_filename.replace('.group', '_') + str(scores) + '.fastq'
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
            if y.isdigit(): # list contains consensus sequence, don't count that one
                grouped_seq += 1
                
    consensusfilename = group_filename.replace('.group', '_consensussequences.fasta') # group consensusfile
    consensusfile = infile.replace('.fasta', '_consensussequences.fasta').replace('.fastq', '_consensussequences.fasta') # total consensusfile
    
    try:  # remove  file if exists
        os.remove(consensusfilename)
    except FileNotFoundError:
        pass 

    MYLOCK.acquire()       
    for j, x in enumerate(grouplist):
        outputfile = group_filename.replace('.group', '_') + str(j) + '.fasta'
        consensusname = group_filename.replace('.group', '_') + str(j) 
        
        t = 0
        for y in x:
            if y.isalpha(): # if it is sequence
                with open(outputfile, 'a') as outputf:
                    outputf.write('>consensus' + '\n' + y + '\n')
                with open(consensusfilename, 'a') as outputf:
                    outputf.write('>consensus_' + str(consensusname) + '(' + str(len(x)-1) + ')\n' + y + '\n')
                with open(consensusfile, 'a') as outputf:
                    outputf.write('>consensus_' + str(consensusname) + '(' + str(len(x)-1) + ')\n' + y + '\n')

            elif y.isdigit():
                t += 1 # count number of sequences in group
        with open('results.tmp', 'a') as rf:
            rf.write('--> ' + outputfile + ' contains ' + str(t) + ' sequences (' + str(round((len(x)-1)*100/len(indexes), 2)) + '%)\n')
    with open('results.tmp', 'a') as rf:
        rf.write(str(grouped_seq) + '/' + str(len(indexes)) + ' sequences assigned in groups for ' 
          + group_filename + ' (' + str(round(grouped_seq*100/len(indexes), 2)) + '%)\n')
    MYLOCK.release()
#==============================================================================
def process_consensuslist(indexes, grouplist, group_filename):  # comparison of sequences with consensus sequences
    global comparelist, len_todolist

    # compare each sequence with consensus1, consensus2,...
    todolist = []
    consensuslist = []
    group_tempfile = group_filename.replace('.group', '.tmp')
#    try:  # remove temporary file if exists
#        os.remove(progressfile)
#    except FileNotFoundError:
#        pass
    
    indexes2 = indexes.copy()  # need a duplicate of indexes to remove items
    for x in grouplist:
        for y in x:
            if y.isdigit(): # list contains consensus sequence, don't count that one
                indexes2.discard(y)  # remove those that are already in a subgroup
                
    for x, y in enumerate(grouplist): # put a number to each consensus that corresponds to the group
        consensuslist.append([x, y[-1]])
        
    comparelist2 = [i for i in comparelist if str(i[3]) in indexes2] # only keep those from group we are working with
#    print('len comparelist2 = ' + str(len(comparelist2)))
#    print('len grouplist = ' + str(len(consensuslist)))
    for x in range(0, len(comparelist2)):
        for y in range(0,len(consensuslist)):
            A1 = comparelist2[x]
            A2 = consensuslist[y]
#            print('A1 = ' + str(A1))
#            print('A2 = ' + str(A2))
            if len(A1[1])*0.7 <= len(A2[1]) <= len(A1[1])*1.3:
                todolist.append([A1,A2])
            else:
                pass  #don't compare if length of sequences differ to much

    len_todolist = len(todolist)
    print(group_filename + '----> ' + str(len(todolist)) + ' comparisons to calculate')
    with open (logfile, 'a') as lf:
        print(group_filename + '--> ' + str(len(todolist)) + ' comparisons to calculate', file = lf)

    try:  # remove temporary file if exists
        os.remove(group_tempfile)
    except FileNotFoundError:
        pass

    templist = []
    q = 0 # number of comparisons done
#    for X in todolist:
    for A1, A2 in todolist:
        q += 1
        try:
            iden = distance(A1[1],A2[1])
            if iden >= similar : # 0.80 is te laag voor 18s
#                    print(str(A1[3]) + '-' + str(A2[0]) + ' : ' + str(iden))
                templist.append((str(A1[3]) + ':' + str(A2[0]) + ':' + str(iden)))
            elif iden < similar:
                iden = distance(A1[1],compl_reverse(A2[1]))
                if iden >= similar:
#                        print(str(A1[3]) + '-' + str(A2[0]) + ':' + str(iden) + ' : ' + ' reverse')
                    templist.append((str(A1[3]) + ':' + str(A2[0]) + ':' + str(iden)))
        except KeyboardInterrupt:
            print("Shutting process down")
#        print(len(templist))
    with open(group_tempfile, 'a') as f:
        for c in templist:
            f.write(c +'\n')
        f.flush()

#==============================================================================        
def update_groups(group_filename, grouplist):
 #    with open (logfile, 'a') as lf:
#        print('Filtering sequences compared with consensus sequences for best hits and add them to the groups', file = lf)
    global comparelist #, grouplist #, consensusgroup
    templist = []
    group_tempfile = group_filename.replace('.group', '.tmp')
    try:
        with open(group_tempfile, 'r') as tf:
            temp = tf.readlines()
        for line in temp:
            templist.append(line.strip().split(':'))
    except FileNotFoundError:
        pass
    templist.sort(key=lambda x: (str(x[0]), str(x[2]))) #sort list based on 1st number (A2) and score 
    for i, j in enumerate(templist[:-1]):
        if j[0] == templist[i+1][0]: # if first index number is the same
            if j[2] <= templist[i+1][2]: # if iden is lower or equal
#                if float(j[2]) < similar:
                j[2] = ''                # mark to remove
    templist = [i for i in templist if i[2] != ''] #only keep those with highest iden
    if len(templist) == 0:
        print(group_filename + '----> no sequences passed criteria')
    else:
#        print('--> Filtering sequences compared with consensus sequences for best hits and add them to the groups')
        print(group_filename + '----> ' + str(len(templist)) + ' sequences added to the groups')

        #..........................................................................
        # Add sequences with high similarity to the groups
        #..........................................................................
        
        with open (logfile, 'a') as lf:
            print('\n*************************************************************', file = lf)
            print('List with similarities > ' + str(similar) , file = lf)
            print('*************************************************************',file = lf)
            for x in templist:
                print(x, file = lf)
                
                
        for x in templist:   
            grouplist[int(x[1])].append(x[0])
    
    
        for i,x in enumerate(grouplist):  # remove consensus from groups
            grouplist[i] = [y for y in x if y.isdigit()]
     
        print(group_filename + '----> Making consensus for each group')
   
        for i,x in enumerate(grouplist):
            consensuslist = []
            for y in x:
                consensuslist.append(comparelist[int(y)][1]) # get the sequence that matches the number
            consensuslist2 = consensus_direction(consensuslist) # get all seq in same direction
            consensus = l_median(consensuslist2)  # create consensuse sequence
            grouplist[i].append(consensus)  # add consensus sequence to the group
    
    try:
        os.remove(group_tempfile)
    except FileNotFoundError:
        pass
    
    return comparelist, grouplist, templist
#==============================================================================
def consensus_direction(consensuslist): #  check if all sequences are in the same direction (F or R)
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
def sort_genes(): # read the input file and sort sequences according to gene
    # remove temporary file if exists
    try:
        filename = infile.replace('.fastq', '_*').replace('.fasta', '_*')
        for x in glob.glob(filename):
            os.remove(x)
        time.sleep(1)
    except FileNotFoundError:
        pass
    read_file(infile)
    process_list(comparelist) 
    update_list(tempfile)
    with open('comparelist.pick', 'wb', buffering=0) as wf:
        pickle.dump(comparelist, wf)

#==============================================================================
def compare_consensus(grouplist):
#    for i,x in enumerate(grouplist):  # only keep consensus from groups
#        grouplist[i] = [y for y in x if not y.isdigit()]
#    for i,[x] in enumerate(grouplist):  # add number to each group
#        grouplist[i] = [i,x]
#    grouplist = [i for i in grouplist if i]  # remove empty sublists    
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
            if len(A1[-1])*1.3 < len(A2[-1]) or len(A2[-1])*1.3 < len(A1[-1]):
                pass  #don't compare if length of sequences differ to much
            else:
                idenlist = []
                iden = distance(A1[-1],A2[-1])
                idenlist.append(iden) # add iden to list
                idenR = distance(A1[-1],compl_reverse(A2[-1]))
                idenlist.append(idenR) # add idenR to list
                idenlist.sort(reverse=True) # sort the list
                iden = idenlist[0] # take the biggest value
                if iden >= 0.96 : 
                    grouplist2[position].extend(A2)
#                    grouplist2[position2] = [] # mark for removal
                elif iden >= 0.90:
                    if len(A1) <= 5 or len(A2) <= 5:
                        grouplist2[position].extend(A2)
#                        grouplist2[position2] = [] # mark for removal
        y += 1  
                 
    grouplist = merge_groups(grouplist2)

    for i,x in enumerate(grouplist):  # remove consensus from groups
        grouplist[i] = [y for y in x if y.isdigit()]
    
    for i,x in enumerate(grouplist):
        consensuslist = []
        for y in x:
            consensuslist.append(comparelist[int(y)][1]) # get the sequence that matches the number
        consensuslist2 = consensus_direction(consensuslist) # get all seq in same direction
        consensus = l_median(consensuslist2)  # create consensuse sequence
        grouplist[i].append(consensus)  # add consensus sequence to the group

##                with open('consensusseq.txt', 'a') as cs:
##                    print('>1\n' + str(A1[-1]) + '\n>2\n' + str(A2[-1]) + '\n' + 
##                          str(iden) +' '+ str(len(A1)) +' '+ str(len(A2)), file = cs)
#            elif iden < 0.90:
#                iden = distance(A1[-1],compl_reverse(A2[-1]))
#                if iden >= 0.90:
#                    with open('consensusseq.txt', 'a') as cs:
#                        print('>1\n' + str(A1[-1]) + '\n>2\n' + str(A2[-1]) + '\n' + 
#                              str(iden) + ' reverse ' + str(len(A1)) +' '+ str(len(A2)), file = cs)
    return grouplist
#==============================================================================   
def rest_reads(indexes, grouplist, group_filename):
    # compare remaining sequences with the consensussequence of the groups,
    # each time with lower similarity, until stable number of sequences in grouplist
    global similar, comparelist
    MYLOCK = Lock()
    print(group_filename + '--> Comparing the rest of the sequences with the consensus sequences')
    print(group_filename + '----> similarity = ' + str(round(similar, 2)))
    process_consensuslist(indexes, grouplist, group_filename)
                
#                  
    comparelist, grouplist, templist = update_groups(group_filename, grouplist)
    min_similar = args.similar_species/100
    while similar >= min_similar: #0.85 is te laag
        while len(templist) > 0:
            process_consensuslist(indexes, grouplist, group_filename)
            comparelist, grouplist, templist = update_groups(group_filename, grouplist)
            grouplist = compare_consensus(grouplist) # compare consensusses with each other 
        else:
            similar = similar - 0.01
            print(group_filename + '----> similarity = ' + str(round(similar, 2)))
            process_consensuslist(indexes, grouplist, group_filename)
            comparelist, grouplist, templist = update_groups(group_filename, grouplist)
            grouplist = compare_consensus(grouplist) # compare consensusses with each other 
            while len(templist) > 0:
                process_consensuslist(indexes, grouplist, group_filename)
                comparelist, grouplist, templist = update_groups(group_filename, grouplist)
          
    grouplist = compare_consensus(grouplist) # compare consensusses with each other 

    MYLOCK.acquire()
    comparelist.sort(key=lambda x: x[3]) #sort list based on index number
    for o, [rec_id, seq, scores, index] in enumerate(comparelist):
        for i, n in enumerate(grouplist):
            if str(index) in n:  # it belongs to a group
                comparelist[o][2] = i      
    MYLOCK.release()
    return grouplist
#==============================================================================
def sort_groups(): # read the gene groups and sort sequences according to species
    global similar, comparelist, resultlist
#    similar = 0.99
    nprocesses = args.nprocesses
    todoqueue = Queue()
    try:
        with open('comparelist.pick', 'rb') as rf:
            comparelist = pickle.load(rf)
    except FileNotFoundError:
        pass 
                   
    for dirpath, dirnames, filenames in os.walk(os.getcwd()):
        for name in filenames:
            if name.endswith('.group'):
                todoqueue.put(name)                    
    try:
        process = [Process(target=sort, args=(todoqueue,)) for x in range(nprocesses)]
        for p in process:
            #ask the processes to stop when all files are handled
            #"STOP" is at the very end of queue
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
    with open('results.tmp', 'r') as rf: # print info of groups at the end
        results = rf.readlines()
        for line in results:
            print(line.strip())
    os.remove('results.tmp')
#==============================================================================
def sort(todoqueue): 
    global similar, comparelist
    consensusfile = infile.replace('.fasta', '_consensussequences.fasta').replace('.fastq', '_consensussequences.fasta') # total consensusfile
    try:  # remove  file if exists
        os.remove(consensusfile)
    except FileNotFoundError:
        pass 
    for group_filename in iter(todoqueue.get, 'STOP'):    #do stuff until infile.get returns "STOP"            
        try:  # remove temporary file if exists
            filename = group_filename.replace('.group', '_*')
            for x in glob.glob(filename):
                os.remove(x)      
    #    sort_genes()
  
            time.sleep(1)
        except FileNotFoundError:
            pass

        similar = 0.95
        indexes, grouplist = read_indexes(group_filename)
        similar = 0.95
        grouplist = rest_reads(indexes, grouplist, group_filename)
        filter_seq(group_filename, grouplist, indexes)

#==============================================================================
"""

- de "unieke" sequenties vergelijken met de groepen-> niet nodig, zijn probleemsequenties (bevatten chimeren en combinaties)

- optie om afzonderlijk in groepen te splitten en nadien afzonderlijk in species te splitten
- oplossing zoeken voor extentie inputfiles (fasta, fastq, fq, fa, ...)
- option om op te slaan in een subfolder
- er zijn nog problemen met het indelen in groepen: sommige genen zitten in verschillende groepen

"""
        
        
        
    
#==============================================================================    
if __name__ == '__main__':
    args = get_arguments()
    infile = args.input
    logfile = infile.replace('.fastq', '.log').replace('.fasta', '.log')
    sort_genes()
    sort_groups()
