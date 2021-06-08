from glob import glob

import gzip
import argparse
import json

import pandas as pd

window = 1000
slide = 10
plot_skew = False

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--folder', help='genome folder', required=True)
args = vars(parser.parse_args())
folder = args['folder']

if __name__ == "__main__":
    files = [f for f in glob(folder + '/*_feature_table.txt.gz')]
    print(files)
    assert len(files) == 1
    file = files[0]

    with gzip.open(file, 'rt') as f:
        df = pd.read_csv(f, sep='\t')

        locs = [(start, end, strand)
                for name, start, end, strand in zip(df.name, df.start, df.end, df.strand)
                if isinstance(name, str) and '16S ribosomal RNA' == name]

        assert(len(locs)) > 0

    with open(file.replace('_feature_table.txt.gz', '_16S.json'), 'w') as fp:
        json.dump(locs, fp)