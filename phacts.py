#!/usr/bin/env python3

import os
import sys
import random
import subprocess
import argparse
from argparse import RawTextHelpFormatter
from subprocess import Popen, PIPE, STDOUT

import numpy as np
from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier

sys.path.pop(0)
import phacts.load_data as load


def unlist(alist):
	return alist[0]

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
	parser.add_argument('-r', '--replicates', type=int, default=10)
	args = parser.parse_args()


	genomes = load.lifestyle()

	labels = dict()
	for genome in genomes.values():
		labels.setdefault(genome.label, []).append(genome.name)
	le = preprocessing.LabelEncoder()
	le.fit( list(labels.keys()) )

	predictions = np.zeros([ args.replicates , len(labels.keys()) ])
	# do ten replicates
	for rep in range(args.replicates):
		# select the genomes
		selected_genomes = list()
		for label in labels:
			selected_genomes.extend(random.sample(labels[label], args.num_genomes))
	
		# select the proteins
		selected_proteins = list()
		for genome in genomes.values():
			for protein in genome.proteins.values():
				if float(protein.importance) > args.cutoff:
					selected_proteins.append(protein)
		selected_proteins = random.sample(selected_proteins, args.num_proteins)
		
		# make the training data
		X = np.zeros([2*args.num_genomes,args.num_proteins])
		y = []
		for i,g in enumerate(selected_genomes):
			genome = genomes[g]
			for j,p in enumerate(selected_proteins):
				X[i,j] = float(genome.similarities.get(p.header, 0))
			y.append(genome.label)
		
		Y = le.transform(y)
		
		clf = RandomForestClassifier(n_estimators=1001)
		clf.fit(X, Y)
	
		prots = ''
		X = np.zeros([1,args.num_proteins])
		for j,p in enumerate(selected_proteins):
			prots += ">temp" + str(j) + "\n"
			prots += p.sequence
			prots += "\n"
			#cmd = "echo '>temp\n" + p.sequence +  "\n' | fasta36 -b 1 -H -q @ " + args.infile + " | grep -m 1 '^Smith-Waterman' | head -n1 | cut -d' ' -f4"
		p = Popen(['fasta35', '-b', '1','@', args.infile], stdout=PIPE, stdin=PIPE, stderr=PIPE)
		output = p.communicate(input= bytes(prots, 'utf-8'))[0]
		flag = False
		j = 0
		for line in output.decode().split('\n'):
			if line.startswith("Library: "):
				flag = True
			elif line.startswith("Smith-Waterman score: ") and flag:
				X[0,j] = float( line.split(' ')[3].replace('%','') )
				j += 1
				flag = False
		
		predictions[rep, :] = unlist(clf.predict_proba(X))
		#label = unlist(le.inverse_transform(np.array([np.argmax(preds)])))
		#predictions.setdefault( label , []).append(preds[np.argmax(preds)])
	print("Class", "probability", "standard deviation", sep ='\t')
	means = predictions.mean(axis=0)
	stdev = predictions.std(axis=0)
	index = np.argmax(means)
	print(unlist(le.inverse_transform([index])), means[index], stdev[index], sep='\t')



