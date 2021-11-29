import lzma
import dill
import os
import pkgutil
import io

def lifestyle():
	data = pkgutil.get_data(__name__, "lifestyle.pkl.xz")
	genomes = None
	f = io.BytesIO(lzma.decompress(data))
	#with lzma.open(data, "rb") as f:
	#	genomes = dill.load(f)
	genomes = dill.load(f)
	return genomes
