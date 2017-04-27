# Python-programming-for-examining-up-and-down-regulated-gene-differences

The code is divided into two parts. The first part takes any file with a column titled “RefSeq Transcript ID” and 
then connects to the NCBI database using the Entrez package. Note that I took out my email address for the entrez system, 
so that should be provided prior to running the code by the user. Here, I provided two files, one with the upregulated oocyte 
genes and one with downregulated genes. For each of the sequence IDs provided, the 3’ UTR are extracted and stored. 
The program then outputs multiple files for analysis and comparison in part 2. The first file includes the length of each 3’ UTR sequence 
and the percentage of As, Ts, Cs, and Gs. The second file gives the number of times each of the 4096 different hexamer motifs occur in the 
sequence. The last file displays the percentage of individual 3’ UTRs in the given file that contain each motif (it does not change based 
on how many occurrence of each motif occur, just whether or not it is present). Lastly, it prints to the screen the percentage of 3’ UTRs 
that had either "AATAAA" and/or "ATTAAA", which are considered the control motifs, and the percentage that had either 
"TTTTAT" and/or "TTTAAT", which would be expected to be larger in the upregulated genes. 


The second python file imports the analysis files for both the up- and down-regulated genes, and uses the pandas package to combine 
these into data frames for comparison. It re-prints from the first section to the screen the percentage of 3’ UTRs that had either 
"AATAAA" and/or "ATTAAA", which are considered the control motifs, and the percentage that had either "TTTTAT" and/or "TTTAAT", which 
would be expected to be larger in the upregulated genes.  Then the statistical differences between the up and down-regulated genes are 
calculated, as well as plots for the various measures. It also displays which motifs are within a certain percentage of the sequences, 
which can be changed by the user. 
