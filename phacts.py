#!/usr/bin/env python3

import os
import sys
import random
import subprocess
import argparse
from subprocess import Popen, PIPE, STDOUT
from sys import platform

import numpy as np
from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier

sys.path.pop(0)
import phacts.load_data as load
#hk.initialize( args.model, os.path.join(params,"parameters_DP09.txt") , os.path.join(params,"multirnafold.conf"), os.path.join(params,"pkenergy.conf") )

def check_fasta(filepath):
	seq = ''
	with open(filepath, mode="r") as f:
		for line in f:
			if not line.startswith(">"):
				seq += line.replace("\n", "").upper()
	ratio = (seq.count('A') + seq.count('C') + seq.count('G') + seq.count('T')) / len(seq)
	return ratio < 0.9

def unlist(alist):
	return alist[0]

def is_valid_file(x):
	if not os.path.exists(x):
		raise argparse.ArgumentTypeError("{0} does not exist".format(x))
	if not check_fasta(x):
		raise argparse.ArgumentTypeError("{0} does not appear to be an amino-acid fasta file".format(x))
	return x

def get_protein_sequences(proteins):
	prots = ''
	for j,p in enumerate(proteins):
		prots += ">temp" + str(j) + "\n"
		prots += p.sequence
		prots += "\n"
	return prots

def parse_arguments():
	usage = '%s [-opt1, [-opt2, ...]] infile' % __file__
	parser = argparse.ArgumentParser(description='', formatter_class=argparse.RawTextHelpFormatter, usage=usage)
	parser.add_argument('infile', type=is_valid_file, help='input file')
	parser.add_argument('-o', '--outfile', action="store", default=sys.stdout, type=argparse.FileType('w'), help='where to write output [stdout]')
	parser.add_argument('-c', '--cutoff', help='Protein importance threshold', type=float, default=0.02)
	parser.add_argument('-g', '--num_genomes', type=int, default=50)
	parser.add_argument('-p', '--num_proteins', type=int, default=600)
	parser.add_argument('-r', '--replicates', type=int, default=10)
	return parser.parse_args()

def spawn_fasta35():
	path = os.path.dirname(load.__file__)
	if platform == "linux" or platform == "linux2":
		path = os.path.join(path, "linux.fasta35")
		return Popen([path, '-b', '1','@', args.infile], stdout=PIPE, stdin=PIPE, stderr=PIPE)
	if platform == "darwin":
		path = os.path.join(path, "osx.fasta35")
		return Popen([path, '-b', '1','@', args.infile], stdout=PIPE, stdin=PIPE, stderr=PIPE)
	
	try:
		p = Popen(['fasta35', '-b', '1','@', args.infile], stdout=PIPE, stdin=PIPE, stderr=PIPE)
		return p
	except:
		raise OSError("known operating system, you will need to manually install fasta35, and make it visible on your PATH")

def select_genomes(labels, args):
	selected_genomes = list()
	for label in labels:
		selected_genomes.extend(random.sample(labels[label], args.num_genomes))
	return selected_genomes

def select_proteins(genomes, args):
	selected_proteins = list()
	for genome in genomes.values():
		for protein in genome.proteins.values():
			if float(protein.importance) > args.cutoff:
				selected_proteins.append(protein)
	selected_proteins = random.sample(selected_proteins, args.num_proteins)
	return selected_proteins

def train_random_forest(genomes, le, selected_genomes, selected_proteins):
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
	return clf

if __name__ == '__main__':
	args = parse_arguments()
	genomes = load.lifestyle()

	labels = dict()
	for genome in genomes.values():
		labels.setdefault(genome.label, []).append(genome.name)
	le = preprocessing.LabelEncoder()
	le.fit( list(labels.keys()) )

	predictions = np.zeros([ args.replicates , len(labels.keys()) ])
	# do ten replicates
	for rep in range(args.replicates):
		selected_genomes = select_genomes(labels, args)
		selected_proteins = select_proteins(genomes, args)
		clf = train_random_forest(genomes, le, selected_genomes, selected_proteins)

		prots = get_protein_sequences(selected_proteins)
		p = spawn_fasta35()
		output = p.communicate(input= bytes(prots, 'utf-8'))[0]
		flag = False
		j = 0
		for line in output.decode().split('\n'):
			if line.startswith("Library: "):
				flag = True
			elif line.startswith("Smith-Waterman score: ") and flag:
				try:
					X[0,j] = float( line.split()[3].replace('%','') )
				except:
					print(output.decode())
					print("The offending line is:")
					print(line)
					exit()
				j += 1
				flag = False
		
		predictions[rep, :] = unlist(clf.predict_proba(X))
		#label = unlist(le.inverse_transform(np.array([np.argmax(preds)])))
		#predictions.setdefault( label , []).append(preds[np.argmax(preds)])
	print(predictions)
	args.outfile.write("Class\tprobability\tstandard deviation\n")
	means = predictions.mean(axis=0)
	stdev = predictions.std(axis=0)
	index = np.argmax(means)
	args.outfile.write("%s\t%s\t%s\n" % (unlist(le.inverse_transform([index])),means[index],stdev[index]) )



