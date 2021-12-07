import lzma
import dill
import os
import pkgutil
import io
import pkg_resources
from subprocess import Popen, PIPE, STDOUT

def lifestyle():
	data = pkgutil.get_data(__name__, "lifestyle.pkl.xz")
	genomes = None
	f = io.BytesIO(lzma.decompress(data))
	#with lzma.open(data, "rb") as f:
	#	genomes = dill.load(f)
	genomes = dill.load(f)
	return genomes

def fasta35():
	p = False
	for f in ['fasta35osx', 'fasta35linux']:
		path = pkg_resources.resource_filename('phacts', f)
		print(path)
		p = Popen([path], stdout=PIPE, stdin=PIPE, stderr=PIPE)
