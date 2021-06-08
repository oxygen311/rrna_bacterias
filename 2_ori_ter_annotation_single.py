from iRep.gc_skew import gc_skew, parse_genomes
from glob import glob

import gzip
import argparse
import json
import os.path

window = 1000
slide = 10
plot_skew = False

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--folder', help='genome folder', required=True)
args = vars(parser.parse_args())
folder = args['folder']

if __name__ == "__main__":
    files = [f for f in glob(folder + '/*_genomic.fna.gz') if not 'rna' in f]
    print(files)
    assert len(files) == 1
    file = files[0]

    if not os.path.isfile(file.replace('_genomic.fna.gz', '_ori_ter.json')):
        with gzip.open(file, 'rt') as fasta_file:
            for name, length, seq in parse_genomes([fasta_file], False):
                ori, ter, _1, _2 = gc_skew(name, length, seq, window, slide, plot_skew)
                # ori, ter = 1, 2
                break  # taking only first contig

        with open(file.replace('_genomic.fna.gz', '_ori_ter.json'), 'w') as fp:
            json.dump([ori, ter, length], fp)