'''
python GetAlleleDist.py -m [MSA/Tree] -i [*.fasta/*.txt] -f 

-m/--mode MSA, tree 

-f/--format fasta, newick 

-i/--input input file 

Purpose: Calculate the genetic distance among tips (particuarly phased alleles) in a multiple sequence alignment or phylogeny. Output the tip(s) with the lowest genetic distance to aid in putative assignment of alleles to subgenomes. 

Requires BioPython 1.81, DendroPy 4.5.2, Numpy 1.20.2, Pandas 1.2.4, phylodm 2.2.1

Jessie Pelosi, University of Florida 2023
'''


import Bio
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
import pandas as pd
import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-m', '--mode', dest = "mode", choices = ['MSA', 'tree'])
parser.add_argument('-i', '--input', dest = "input")
parser.add_argument('-f', '--format', dest = "format")
args = parser.parse_args()

MODE = args.mode
INPUT = args.input
FORMAT = args.format

if MODE == "MSA":
	# Read alignment
	aln = AlignIO.read(open(INPUT), FORMAT)

	# Read record names into empty object 
	recs = []
	for record in aln:
	     recs.append(record.id)

	# Calculate genetic distance between all records in alignment 
	calculator = DistanceCalculator('identity')
	dm = calculator.get_distance(aln)

	# Convert genetic distance matrix to dataframe, add in names of records 
	df = pd.DataFrame(list(dm))
	df.columns = recs
	df.index = recs

	# Replace diagonals with NaN, this will prevent self-self zero distances 
	np.fill_diagonal(df.values, np.NaN)

	# Write allele distances to csv file 
	df.to_csv(INPUT + ".DistMatrix.csv")

	# Get the lowest genetic distances for each record, allow for equal distances 
	s = df.stack()
	bestAllele = (s[s.eq(s.groupby(level=0).transform('min'))]
	  .reset_index()
	  .groupby('level_0')['level_1'].apply(list)
	 )

	# Convert to dataframe and write to csv 
	BAdf = pd.DataFrame(bestAllele)
	BAdf.columns = ['BestHits']
	BAdf.to_csv(INPUT + ".ClosestAllele.csv")

elif MODE == "tree":
	from phylodm import PhyloDM

	# Import tree from a DendroPy
	import dendropy
	tree = dendropy.Tree.get(
	    path=INPUT,
	    schema=FORMAT)
	pdm = PhyloDM.load_from_dendropy(tree)
	# Calculate the PDM
	dm = pdm.dm(norm=False)
	labels = pdm.taxa()

	# Convert distance matrix to dataframe 
	df = pd.DataFrame(list(dm))
	df.columns = labels
	df.index = labels

	# Replace diagonals with NaN, this will prevent self-self zero distances 
	np.fill_diagonal(df.values, np.NaN)

	# Write allele distances based on tree to csv file 
	df.to_csv(INPUT + "DistMatrix.csv")

	# Get the lowest genetic distances for each record, allow for equal distances 
	s = df.stack()
	bestAllele = (s[s.eq(s.groupby(level=0).transform('min'))]
	  .reset_index()
	  .groupby('level_0')['level_1'].apply(list)
	 )

	# Convert to dataframe and write to csv 
	BAdf = pd.DataFrame(bestAllele)
	BAdf.columns = ['BestHits']
	BAdf.to_csv(INPUT + ".ClosestAlleles.csv")

else:
	print("Missing mode (-m) argument!")

