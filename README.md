# amplicon_sorter

Sorts amplicons from Nanopore sequencing data based on similarity, creates groups with highly similar sequences and create the consensus sequence for each group of sequences.

(Only works on Linux because it uses multiprocessing)

(This is a work in progress, but already producing good enough results to share the script.)

Requirements:

-   Python 3
-   python3-dev, python3-setuptools, python3-pip, python3-wheel (`sudo apt-get install python3-dev python3-setuptools python3-pip python3-wheel`)
-   c-implementation of Levenshtein: [https://pypi.org/project/python-Levenshtein/](https://pypi.org/project/python-Levenshtein/) (`python3 -m pip install python-Levenshtein`)
-   biopyton (`sudo apt-get install python3-biopython`)
-   matplotlib (`sudo apt-get install python3-matplotlib`)

### Options:

`-i, --input`: Input file in fastq or fasta format (auto-detect). Make sure the inputfile is named as .fasta or .fastq because it replaces the extension in parts of the script.

`-min, --minlength`: Minimum readlenght to process. Default=300

`-max', '--maxlength`: Maximum readlenght to process. Default=No limit

`-maxr', '--maxreads`: Maximum number of reads to process. Default=10000

`-np', '--nprocesses`: Number of processors to use. Default=1

`-sg', '--similar_genes`: 'Similarity to sort genes in groups (value between 50 and 100). Default=55.0

`-ssg', '--similar_species_groups`: Similarity to CREATE species groups (value between 50 and 100). Default=92.0

`-ss', '--similar_species`: Similarity to ADD sequences to a species group (value between 50 and 100). Default=85.0

`-sfq', '--save_fastq`: Save the results also in fastq files (fastq files will not contain the consensus sequence)

`-ra', '--random`: Takes random reads from the inputfile if --maxreads is lower than total number of reads that passed criteria

`-of', '--outputfolder`: Save the results in the specified outputfolder. Default = current working directory

`-ho', '--histogram_only`: Only makes a read length histogram.  Can be interesting to see what the minlength and maxlength setting should be.

`-so', '--species_only`: Only creates species groups and sort to species level.  This can only be done if the whole script has run once without this option.  This is to save time if you play with the `--similar_species` and/or `--similar_species_groups` parameters.  It is using the same `--maxreads`, `--minlength`, `--maxlength` data that is produced in the first part of the script, so those 3 parameters are ignored here.

### How it works (in short):

1.  The script reads the inputfile and creates a read length histogram of all the reads in the file.
2.  It starts processing a selection of the reads (based on minlength, maxlenght, maxreads). It saves the result in .group files that contain reads of the same gene (e.g. group_1 with 18S reads, group_2 with COI reads, group_3 with ITS reads).
3.  It processes the group files to sort out the genes to species or genus level and saves this to different files. (e.g. file_1_1.fasta is 18S from species1, file1_2.fasta is 18S from species2, ...)
4.  Each outputfile contains at the end the consensus file of the sequences in the file.
5.  Files are produced per group with all consensus sequences for that group. A file is produced with all consensus sequences of all groups.
6.  Files are produced with 'unique' sequences (script does not find a group where it belongs to according to the settings)


![amplicon_sorter](https://github.com/avierstr/amplicon_sorter/blob/master/ampliconsorter.jpg)

### Command examples:

`python3 amplicon_sorter.py -i inputset1.fastq` : process file with default settings.

`python3 amplicon_sorter.py -i inputset1.fastq --similar_species 80 --nprocesses 8`: process file on 8 cores with allowing the similarity between sequences of the same gene as low as 80%

### My experience so far:

-   Guppy produces reads with Q-score >=7 by default. If I filter reads to Q-score >= 10 or Q-score >= 11, the end result (similarity % of the consensus to a species in BLAST) is for all comparable. The only difference is that I get more reads per species group for Q-score 11 than for Q-score 7, and that there are less reads in the 'unique' groups.
-   The same results if I lower the `--similar_species` value from 90% to 80%. Above 90% the similarity in Blast is getting lower, between 80 and 90 the similarity is comparable. The lower I go, the more sequences in the species groups.
-   10000 reads (`--maxreads`) on 12 cores takes several hours to finish.

### Todo:

-  For now, only files with .fasta or .fastq extensions are allowed. Make it foolproof to allow other extension as long it contains sequences in fastq or fasta format.
-  Try to improve speed for comparison of sequences.

### Release notes:

2020/4/27:

- Little change in step 3 that improves the processing speed.
- Made names of temporary files unique in case you want to process different files in one folder.

2020/4/17:

- Added vertical lines to the read length histogram with the minimum and maximum read length that were given as read length selections.![readlengthhistogram](https://github.com/avierstr/amplicon_sorter/blob/master/readlengthhistogram.jpg)
- Option to create a read length histogram only (`--histogram_only`).
- Option to change the similarity when CREATING species groups (`--similar_species_groups`).
- Option (`--species_only`) to play with the `--similar_species` and `--similar_species_groups` parameters without having to start all over.  

2020/4/3:

 - Option to save results in a subfolder (`--outputfolder`).
 - More little improvement to sort out the groups in step 2. (The result is still not what I had in mind: sometimes reads of ITS and 18S are in the same group.  The reason for this is that ITS can be very different from species to species.  If I set similarity to 55%, it can put some ITS in the 18S group.  If I set it higher, than COI is not in one group.  This has no influence on the end results in step 3-4). 

2020/3/12:

-   Option to take a random selection of the reads from the inputfile (`--random`).
-   Little improvement to sort out the groups in step 2.

2020/3/5:

-   Fasta or fastq files possible as input (autodetect).
-   Fastq files as output option (when input is fastq) (`--save_fastq`).

