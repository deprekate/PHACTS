import lzma
import dill
import pkgutil
import io
import pkg_resources

def lifestyle():
	data = pkgutil.get_data(__name__, "lifestyle.pkl.xz")
	genomes = None
	f = io.BytesIO(lzma.decompress(data))
	genomes = dill.load(f)
	return genomes