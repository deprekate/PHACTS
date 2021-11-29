import os
import sys
import random
import pickle
import dill
import lzma


class Protein(dict):
	def __init__(self, header=None, sequence=''):
		self.header = header
		self.sequence = sequence
		self.importance = 0

class Genome(dict):
	def __init__(self, name=None, label=None):
		self.name = name
		self.label = label
		self.proteins = dict()
		self.similarities = dict()


# load genomes
genomes = dict()
fp = open("classes_lifestyle")
for line in fp:
	name,label = line.rstrip().split('\t')
	genome = Genome(name, label)
	protein = Protein()
	fp = open("gen/" + name)
	for line in fp:
		if line.startswith(">"):
			genome.proteins[protein.header] = protein
			protein = Protein(header=line[1:].split()[0])
		else:
			protein.sequence += line.rstrip()
	genome.proteins[protein.header] = protein
	genomes[name] = genome
fp.close()

#load importances
fp = open("classes_lifestyle.importance")
for line in fp:
	p,value = line.rstrip().split('\t')
	pg = p.split('.peg')[0][4:] + '.fas'
	genomes[pg].proteins[p].importance = value
fp.close()

fp = open("sims.txt")
for line in fp:
	key,value = line.rstrip().split('\t')
	p,g = key.split('<->')
	pg = p.split('.peg')[0][4:] + '.fas'
	if pg in genomes and g in genomes:
		genomes[g].similarities[p] = value
fp.close()

print(genomes["101570.3.fas"].similarities["fig|101570.3.peg.1"])

with lzma.open("lifestyle.pkl.xz", "wb") as f:
    #pickle.dump(genomes, f, pickle.DEFAULT_PROTOCOL)
    dill.dump(genomes, f, pickle.DEFAULT_PROTOCOL)
'''
data = None
with lzma.open("lifestyle.pkl.xz", "rb") as f:
	data = pickle.load(f)

print("BEGIN")
selected = list()
keys = list(data.keys())
while( len(selected) < 10):
	key = random.choice(keys)
	if float(data[key].importance) > 0.001:
		selected.append(key)
print(selected)
'''
