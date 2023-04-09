'''
usage: python distribute_alleles.py -i ORIG.fasta -p PHASED.fasta -s SUBGENOMES.csv -l LOCUS

-l: locus name
-i: Original fastas with phased taxa removed 
-p: Phased fastas with phased taxa included 
-s: CSV with subgenome assignments 
	Subgenomes database should contain three columns: 
		1) the locus which matches the name of the locus used for the fastas
		2) the sample name appended with __REF, __1, or __2 (deafult output from PATE) 
		3) the assigned subgenome (A, B). If left blank, the sample will not be included in the new fasta with phased alleles assigned to subgenomes. 

purpose: distribute alleles that are assigned to subgenomes into the appropriate locus and rename with subgenome assignment

'''

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest = "input")
parser.add_argument('-p', '--phased', dest = "phased")
parser.add_argument('-s', '--subgenomes', dest = "subgenomes")
parser.add_argument('-l', '--locus', dest = "locus")
args = parser.parse_args()

FASTA = args.input
PHASED = args.phased
SUBGENOME = args.subgenomes
LOCUS = args.locus

import pandas as pd
import re

# Read in subgenomes dataframe 
subs = pd.read_csv(SUBGENOME, header = 0)

# Subset the dataframe to only include the locus of interest 
locus_df = subs[subs.locus == LOCUS].dropna()

# Iterate through the new dataframe sample by sample, and replace the 
from Bio import SeqIO
for record in SeqIO.parse(PHASED, "fasta"):
    for index, row in locus_df.iterrows():
        sample = row[0]
        sub = row[2]
        if record.name in sample:
            record.name = re.sub(r"__[A-Z0-9]*", "_"+str(sub), record.name)
            with open(FASTA, 'a') as fasta:
                fasta.write("\n>" + record.name + '\n' + str(record.seq) + "\n")
