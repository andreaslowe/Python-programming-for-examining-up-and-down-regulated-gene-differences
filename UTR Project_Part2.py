#This code reads in the analysis files created in the UTR Project Part 1 code and performs statistics and graphical comparisons
#between the up- and down-regulated genes

#import analysis files
import pandas as pd  
MotifCounts_Down = pd.read_csv("MotifCounts_Downregulated oocyte genes.csv", header = None)
MotifPercent_Down = pd.read_csv("MotifPercents_Downregulated oocyte genes.csv", header = None)
Nucleotides_Down = pd.read_csv("Nucleotides_Downregulated oocyte genes.csv")
MotifCounts_Up = pd.read_csv("MotifCounts_Upregulated oocyte genes.csv", header = None)
MotifPercent_Up = pd.read_csv("MotifPercents_Upregulated oocyte genes.csv", header = None)
Nucleotides_Up = pd.read_csv("Nucleotides_Upregulated oocyte genes.csv")

#check for control and upregulated motifs, sum number of UTRs that contain these motifs and determine percentage
control_motif_down = 0
upreg_control_motif_down = 0

control_motif_up = 0
upreg_control_motif_up = 0

MotifCounts_Down.set_index(0)
control_down = MotifCounts_Down.loc[(MotifCounts_Down[0] == "AATAAA") | (MotifCounts_Down[0] == "ATTAAA")]
control_down = control_down.drop(0, axis = 1)
for column in control_down:
    if sum(control_down[column]) >= 1:
            control_motif_down += 1
control_percent_down = 100*(control_motif_down/float(len(control_down.columns)))
            
upreg_control_down = MotifCounts_Down.loc[(MotifCounts_Down[0] == "TTTTAT") | (MotifCounts_Down[0] == "TTTAAT")]
upreg_control_down = upreg_control_down.drop(0, axis = 1)
for column in upreg_control_down:
    if sum(upreg_control_down[column]) >= 1:
            upreg_control_motif_down += 1
upreg_percent_down = 100*(upreg_control_motif_down/float(len(control_down.columns)))
            
MotifCounts_Up.set_index(0)
control_up = MotifCounts_Up.loc[(MotifCounts_Up[0] == "AATAAA") | (MotifCounts_Up[0] == "ATTAAA")]
control_up = control_up.drop(0, axis = 1)
for column in control_up:
    if sum(control_up[column]) >= 1:
            control_motif_up += 1
control_percent_up = 100*(control_motif_up/float(len(control_up.columns)))
            
upreg_control_up = MotifCounts_Up.loc[(MotifCounts_Up[0] == "TTTTAT") | (MotifCounts_Up[0] == "TTTAAT")]
upreg_control_up = upreg_control_up.drop(0, axis = 1)
for column in upreg_control_up:
    if sum(upreg_control_up[column]) >= 1:
            upreg_control_motif_up += 1
upreg_percent_up = 100*(upreg_control_motif_up/float(len(control_up.columns)))
            
#prepare to plot sequence lengths and percentage of A, T, C, and G
seq_len = pd.DataFrame.from_items([("Up Regulated", Nucleotides_Up["Sequence Length"]), ("Down Regulated", Nucleotides_Down["Sequence Length"])])
A_percent = pd.DataFrame.from_items([("Up Regulated", Nucleotides_Up["A Percent"]), ("Down Regulated", Nucleotides_Down["A Percent"])])
T_percent = pd.DataFrame.from_items([("Up Regulated", Nucleotides_Up["T Percent"]), ("Down Regulated", Nucleotides_Down["T Percent"])])
C_percent = pd.DataFrame.from_items([("Up Regulated", Nucleotides_Up["C Percent"]), ("Down Regulated", Nucleotides_Down["C Percent"])])
G_percent = pd.DataFrame.from_items([("Up Regulated", Nucleotides_Up["G Percent"]), ("Down Regulated", Nucleotides_Down["G Percent"])])

#Code to plot parameters - copy into command line to see
seq_len.plot(kind = "box")
A_percent.plot(kind = "box") 
T_percent.plot(kind = "box")
C_percent.plot(kind = "box")
G_percent.plot(kind = "box")
#**figures must be explicitly closed using: matplotlib.pyplot.close("all")

#create dataframe that compares motifs in percent of sequences Up vs Down
MotifPercent = pd.DataFrame.from_items([("Motifs", MotifPercent_Up.ix[:,0]), ("Up Regulated", MotifPercent_Up.ix[:,1]), ("Down Regulated", MotifPercent_Down.ix[:,1])])

#stats differences in motif counts
import scipy.stats as stats
counts_down = MotifCounts_Down.drop(0, axis=1)
counts_down = counts_down.T
counts_up = MotifCounts_Up.drop(0, axis=1)
counts_up = counts_up.T

p_value = []
index = []
for columns in counts_down:
    t, p = stats.ttest_ind(counts_down[columns], counts_up[columns], nan_policy= "omit")
    if p <= 0.05:
        index.append(columns)
        p_value.append(p)
significant_motifs = MotifCounts_Down.iloc[index][0]

#Print percentage of control and upregulated
print "Percent control motifs (AATAAA or ATTAAA) in down-regulated genes: %f" %control_percent_down
print "Percent upregulated motifs (TTTTAT or TTTAAT) in down-regulated genes: %f" % upreg_percent_down
print "Percent control motifs (AATAAA or ATTAAA) in up-regulated genes: %f" %control_percent_up
print "Percent upregulated motifs (TTTTAT or TTTAAT) in up-regulated genes: %f" % upreg_percent_up

#To look at different cutoffs for the motifs that are in greater than 60% of sequences (can change number)
print "Motifs in greater than 60 percent of up regulated:" 
print MotifPercent_Up[MotifPercent_Up[1] > 55]
print "Motifs in greater than 60 percent of down regulated:"
print MotifPercent_Down[MotifPercent_Down[1] > 55]

#standard t-test for sequence length, can repeat for other factors - should probably look into chi-squared test and correcting for multiple comparisons (Type I error), etc
print "statistical differences in sequence length:"
print stats.ttest_ind(seq_len["Up Regulated"], seq_len["Down Regulated"], nan_policy= "omit")

print "statistical differences in A_percent:"
print stats.ttest_ind(A_percent["Up Regulated"], A_percent["Down Regulated"], nan_policy= "omit")

print "statistical differences in T_percent:"
print stats.ttest_ind(T_percent["Up Regulated"], T_percent["Down Regulated"], nan_policy= "omit")

print "statistical differences in C_percent:"
print stats.ttest_ind(C_percent["Up Regulated"], C_percent["Down Regulated"], nan_policy= "omit")

print "statistical differences in G_percent:"
print stats.ttest_ind(G_percent["Up Regulated"], G_percent["Down Regulated"], nan_policy= "omit")

print "There were %i motifs that had significant differencs in the counts between up and down-regulated genes." %len(index) 
print "Type significant_motifs to see them and p_value to see associated values"
