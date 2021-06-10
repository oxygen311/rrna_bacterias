from iRep.gc_skew import gc_skew, parse_genomes
from glob import glob

import gzip
import argparse
import json
import os.path

import numpy as np

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
            parsed = list(parse_genomes([fasta_file], False))
            max_len = max(len for _1, len, _2 in parsed)

            for name, length, seq in parsed:
                if length == max_len:
                    ori, ter, skew, c_skew = gc_skew(name, length, seq, window, slide, plot_skew)
                    ori_i, ter_i = ori // slide, ter // slide

                    if ori < ter:
                        repl_1 = np.mean(skew[1][ori_i:ter_i])
                        repl_2 = np.mean(skew[1][ter_i:] + skew[1][:ori_i])
                    else:
                        repl_2 = np.mean(skew[1][ter_i:ori_i])
                        repl_1 = np.mean(skew[1][ori_i:] + skew[1][:ter_i])

                    break  # taking only first contig

        with open(file.replace('_genomic.fna.gz', '_ori_ter.json'), 'w') as fp:
            json.dump([name, ori, ter, repl_1, repl_2, np.min(c_skew[1]), np.max(c_skew[1]), length], fp)