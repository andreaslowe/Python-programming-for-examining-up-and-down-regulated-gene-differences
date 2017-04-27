#This code analyzes the up and down-regulated oocyte genes, by extracting the sequence from the NCBI database and then couting motifs and other markers 
#from the 3' UTR region 
 
file_id = "Upregulated oocyte genes.csv"  #also run for the "Downregulated oocyte genes.csv" file

sequence_ids = []
import pandas as pd
data = pd.read_csv(file_id)
data = data['RefSeq Transcript ID'].dropna().tolist()
for i in range(len(data)):
    sequence_ids.append(data[i].split("///")[0])  #takes the first refseqid in the list (since the rest are different names for the same gene)

seq_num = len(sequence_ids)
UTR_list = []
id_list = []

#connect to NCBI database and extract coding region
from Bio import Entrez, SeqIO
Entrez.email = "" #**** should always tell the NCBI who you are - I took out my email address so that anyone can use**** 
handle = Entrez.efetch(db="nuccore", id=sequence_ids, rettype="gb", retmode="text")
records = SeqIO.parse(handle,"genbank")
for record in records:
    sequence = "%s" % record.seq
    id_list.append(record.id) 
    for feature in record.features:
        if feature.type =="CDS":
            UTR_loc = feature.location.end+1 #find location of the end of the coding region
            UTR = sequence[UTR_loc:] #extract 3'UTR
            #only store if the UTR has at least 1 nucleotide
            if len(UTR) > 0: 
                UTR_list.append(UTR) #append to list
handle.close()

#store number of sequences
seq_num = len(UTR_list)

#create list of motif dictionaries and empty values for lists used
dict_motifs = [dict() for x in range(seq_num)]
seq_lengths = []
A_percent = []
T_percent = []
G_percent = []
C_percent = []
max_count = [] 
max_key = [0] * seq_num
control_motif = 0
upreg_control_motif = 0

#add motifs counts for each of the 3'UTRs to the empty dictionaries, store the length and count of each individual nucleotide, add number of sequences with control motifs
for x in range(seq_num):    
    
    seq_lengths.append(len(UTR_list[x]))
    A_percent.append(100*UTR_list[x].count("A")/float(len(UTR_list[x]))) 
    T_percent.append(100*UTR_list[x].count("T")/float(len(UTR_list[x]))) 
    G_percent.append(100*UTR_list[x].count("G")/float(len(UTR_list[x])))
    C_percent.append(100*UTR_list[x].count("C")/float(len(UTR_list[x])))     
    
    #enumerate hexamer motifs                         
    for letter1 in ["A", "T", "G", "C"]:
        for letter2 in ["A", "T", "G", "C"]:
            for letter3 in ["A", "T", "G", "C"]:
                for letter4 in ["A", "T", "G", "C"]:
                    for letter5 in ["A", "T", "G", "C"]:
                        for letter6 in ["A", "T", "G", "C"]:
                            motif = letter1 + letter2 + letter3 + letter4 + letter5 + letter6 
                            dict_motifs[x][motif] = UTR_list[x].count(motif)

    #return list of most frequent motif(s) in each sequence
    max_count.append(max(dict_motifs[x].values()))
    max_key[x] = [k for k, v in dict_motifs[x].items() if v == max_count[x]]             
                
    #check for control and upregulated motifs, sum number of UTRs that contain these motifs
    if ((dict_motifs[x]["AATAAA"] >= 1) | (dict_motifs[x]["ATTAAA"] >= 1)):  #should be in most 
        control_motif += 1
    if ((dict_motifs[x][("TTTTAT")] >= 1) | (dict_motifs[x]["TTTAAT"] >= 1)): #should be in upregulated
        upreg_control_motif += 1

#calculate the percentage of sequences that contain the control and upregulated motif   
control_percent = 100*(control_motif/float(seq_num))
upreg_percent = 100*(upreg_control_motif/float(seq_num))

#convert to data frame
motif_counts = pd.DataFrame(dict_motifs)
motif_counts = motif_counts.T #transpose

#save counts to csv file
filename = "MotifCounts_" + file_id
motif_counts.to_csv(filename, header = False)

#make list of the percentage of sequences that contain each motif, and second output of only present motifs
motif_tmp = motif_counts
motif_tmp[motif_tmp>1] = 1
motif_percents = 100 * (motif_tmp.sum(axis = 1)/seq_num)
present_motifs = motif_percents[motif_percents > 0]

#create data frames of the percentage of motifs and the nucleotide makeup
percent_motifs = pd.DataFrame(present_motifs)
nucleotides = pd.DataFrame.from_items([('Sequence Length', seq_lengths), ('A Percent',A_percent),('T Percent',T_percent),('C Percent',C_percent),('G Percent',G_percent)])

#save other descriptors to file
filename2 = "MotifPercents_" + file_id
filename3 = "Nucleotides_" + file_id
percent_motifs.to_csv(filename2, header = False)
nucleotides.to_csv(filename3, index = False)

#consider using present_motifs seq_lengths, A_percent (etc), max_count, max_key for comparisons 
print "For the file %s there are: " %file_id
print "Percent control motifs (AATAAA or ATTAAA): %f" %control_percent
print "Percent upregulated motifs (TTTTAT or TTTAAT): %f" % upreg_percent
print "Number of motifs in all sequences: %d" %len(present_motifs)
