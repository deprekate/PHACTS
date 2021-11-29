#!/usr/bin/env python3

import os
import sys
import random
import pickle
import lzma
import subprocess
import argparse
from argparse import RawTextHelpFormatter

import dill
import numpy as np
from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier


sys.path.pop(0)
import phacts.load_data as load

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	return x

if __name__ == '__main__':
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', type=is_valid_file, help='input file')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write output [stdout]')
	parser.add_argument('-c', '--cutoff', help='The minimum cutoff length for runs', type=int, default=0)
	parser.add_argument('-g', '--num_genomes', type=int, default=50)
	parser.add_argument('-p', '--num_proteins', type=int, default=600)
	args = parser.parse_args()

	#genomes = None
	#with lzma.open("lifestyle.pkl.xz", "rb") as f:
	#	genomes = dill.load(f)

	genomes = load.lifestyle()
	
	# select the genomes
	labels = dict()
	selected_genomes = list()
	for genome in genomes.values():
		labels.setdefault(genome.label, []).append(genome.name)
	
	for label in labels:
		selected_genomes.extend(random.sample(labels[label], args.num_genomes))
	
	# select the proteins
	selected_proteins = list()
	for genome in genomes.values():
		for protein in genome.proteins.values():
			if float(protein.importance) > args.cutoff:
				selected_proteins.append(protein)
	selected_proteins = random.sample(selected_proteins, args.num_proteins)
	
	
	X = np.zeros([2*args.num_genomes,args.num_proteins])
	y = []
	for i,g in enumerate(selected_genomes):
		genome = genomes[g]
		for j,p in enumerate(selected_proteins):
			X[i,j] = float(genome.similarities.get(p.header, 0))
		y.append(genome.label)
	
	le = preprocessing.LabelEncoder()
	le.fit(y)
	Y = le.transform(y)
	
	clf = RandomForestClassifier(max_depth=2, random_state=0)
	clf.fit(X, Y)
	
	X = np.zeros([1,args.num_proteins])
	for j,p in enumerate(selected_proteins):
		cmd = "echo '>temp\n" + p.sequence +  "\n' | fasta36 -b 1 -H -q @ " + args.infile + " | grep '^Smith-Waterman' | head -n1 | cut -d' ' -f4"
		pid = subprocess.getoutput(cmd)[:-1]
		X[0,j] = float(pid)
	
	preds = clf.predict_proba(X)
	label = le.inverse_transform(np.array([np.argmax(preds)]))
	print(label, preds)




