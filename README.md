
# amplicon_sorter

A tool for reference-free sorting of ONT sequenced amplicons based on their similarity in sequence and length and for building solid consensus sequences.
The limit for separating closely related species within a sample is currently around 95 - 96%.

For more detailed explanation, please read [Amplicon_sorter_manual.pdf](https://github.com/avierstr/amplicon_sorter/blob/master/Amplicon_sorter_manual.pdf).

(Only works on Linux because it uses multiprocessing)

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

`-sc', '--similar_consensus`: Similarity to COMBINE groups based on the consensus sequence (value between 50 and 100). Default=96.0

`-sfq', '--save_fastq`: Save the results also in fastq files (fastq files will not contain the consensus sequence)

`-ra', '--random`: Takes random reads from the inputfile.  The script does NOT compare al sequences with each other, it compares batches of 1.000 with each other.  You can use this option and sample reads several times and compare them with other reads in other batches.  So it is possible to have an inputfile with 10.000 reads and sample random 20.000 reads from that inputfile.  The script will run 20 batches of 1.000 reads.  This way, the chance to find more reads with high similarity is increasing when there are a lot of different amplicons in the sample.  No need to do that with samples with 1 or 2 amplicons.

`-a', '--all`: Compare all selected reads with each other.  Only advised for a small number of reads (< 10000)

`-o', '--outputfolder`: Save the results in the specified outputfolder. Default = current working directory

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
### Workflow:

Filter your inputfile for reads >= Q12 with NanoFilt ([https://github.com/wdecoster/nanofilt](https://github.com/wdecoster/nanofilt)) or other quality filtering software.  Use that Q12 inputfile for Amplicon_sorter. (Lower quality reads can be used but will result in longer processing time and a lower percentage of reads that will assigned to a species.

Copy the Amplicon_sorter.py script in the same folder as your inputfile.

### Command examples:
*Produce a read length histogram of your inputfile:*

`python3 amplicon_sorter.py -i infile.fastq –o outputfolder -min 650 -max 750 -ho`: produce the readlength histogram of infile.fastq in folder outputfolder.  This gives you the information on the number of reads between 650 and 750 bp.

*Sample with one species amplicon of 750 bp:*

`python3 amplicon_sorter.py -i infile.fastq -o outputfolder -np 8 -min 700 -max 800 -maxr 1000`

process infile.fastq with default settings, save in folder outputfolder, run on 8 cores, minimum length of reads = 700, max length of reads = 800, use 1000 reads.  This will sample the first 1000 reads between 700 and 800 bp of the inputfile.  If you add the -ra (random) option to the command line, it will sample 1000 random reads between 700 and 800 bp.

*Sample with 2 species: an amplicon of 700 bp and one of 1200 bp:*

`python3 amplicon_sorter.py -i infile.fastq -o outputfolder -np 8 -min 650 -max 1250 -maxr 2000`  

*Metagenetic sample with several amplicons between 600 and 3000 bp, unknown number of species, 30000 reads in the inputfile:*

`python3 amplicon_sorter.py -i infile.fastq -o outputfolder -np 8 -min 550 -max 3050 -maxr 30000`

*Metagenetic sample with several amplicons between 600 and 3000 bp, unknown number of species, 30000 reads in the inputfile, one low abundant species (< 2% reads):*

`python3 amplicon_sorter.py -i infile.fastq -o outputfolder -np 8 -min 550 -max 3050 -ra -maxr 600000`

By random sampling 20x the maximum number of reads, it is possible to find low abundant species.

### Todo:

- Try to improve speed for comparison of sequences.
- Try to fine tune sorting to species level.  If there are species in the sample that are more than 94% similar, the script has difficulties to separate them.  Species with a higher similarity are often grouped together and give a lower consensus sequence because it is the "average" of 2 or more closely related species.
- Improve random sampling to sequences of the same length and not to all reads.
- Check for bug: there is sometimes an error in the percentage of reads assigned per group (sometimes > 100%).  This has no effect on the sorting or consensus made, only on the information how many reads are assigned.

### Release notes:
2021/09/21:
-fixed bug that is important for sorting closely related species. (was wrongfully removed in previous version)

2021/09/11:
- speed and memory improvement when sorting the compared reads for best matches
- added option `--similar_consensus` 
- added option `--all`
- save all parameters in results.txt to know afterwards with which parameters a run is performed.

2021/08/08:
- speed improvement by changing allowed difference in length from 1.1 to 1.05 (less comparisons to be done) and 1.08 for comparison of consensus sequences
- change default of `--similar_species_groups` from 0.92 to 0.93

2021/07/13:
- speed improvement by changing allowed difference in length from 1.3 to 1.1 (less comparisons to be done)
- little changes that improve speed a bit

2021/05/28:
- improvement on combining groups
- small speed improvement when creating consensus to compare groups

2021/05/19:
- small speed improvement when comparing remaining sequences with consensus in groups

2021/05/13:
- chunck number of comparisons in groups of 500K to save memory and process parallel
- minor bug fixes and small speed improvements
- improvements in number of groups created

2021/05/06:
- If no max number of reads is given, use all reads
- write version of script in results.txt
- write number of reads in length selection to readlength histogram
- combine groups if consensus is at least 99% similar
- save time by making consensus of max 500 random reads in group
- major speed improvement when comparing a lot of sequences.  When the number of sequences is x 100, calculation time is x 10.000.  Subsampling in batches of 1000 sequences decreases calculation time x 10 (it does NOT compare al sequences with each other, it compares batches of 1000 with each other) 

2020/5/20:
- Save a file "results.txt" in outputfolder with how many reads are in which file.
- Catching error when there are less reads available than asked in `--maxreads`.

2020/5/6:

- Little cosmetic change in read length histogram.
- Corrected a minor bug with (`-ho --histogram_only`) option.

2020/4/27:

- Little change in step 3 that improves the processing speed.
- Made names of temporary files unique in case you want to process different files in one folder.

2020/4/17:

- Added vertical lines to the read length histogram with the minimum and maximum read length that were given as read length selections.![readlengthhistogram](https://github.com/avierstr/amplicon_sorter/blob/master/readlengthhistogram.jpg)
- Option to create a read length histogram only (`--histogram_only`).
- Option to change the similarity when CREATING species groups (`--similar_species_groups`).
- Option (`--species_only`) to play with the `--similar_species` and `--similar_species_groups` parameters without having to start all over.  

2020/4/3:

 - Option to save results in a subfolder (`-o --outputfolder`).
 - More little improvement to sort out the groups in step 2. (The result is still not what I had in mind: sometimes reads of ITS and 18S are in the same group.  The reason for this is that ITS can be very different from species to species.  If I set similarity to 55%, it can put some ITS in the 18S group.  If I set it higher, than COI is not in one group.  This has no influence on the end results in step 3-4). 

2020/3/12:

-   Option to take a random selection of the reads from the inputfile (`--random`).
-   Little improvement to sort out the groups in step 2.

2020/3/5:

-   Fasta or fastq files possible as input (autodetect).
-   Fastq files as output option (when input is fastq) (`--save_fastq`).





