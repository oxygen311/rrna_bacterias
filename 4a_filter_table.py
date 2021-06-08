import pandas as pd
import json

table = 'all_refseq_bacteria.tsv'


def process_row_ori_ter(accession, asm, refseq_category, organism_name, infraspecific_name):
    ori_ter_file = f'refseq/bacteria/{accession}/{accession}_{asm}_ori_ter.json'
    with open(ori_ter_file) as f:
        data = json.load(f)
        print(data)


if __name__ == "__main__":
    ncbi_down_df = pd.read_csv(table, sep='\t')

    ncbi_down_df = ncbi_down_df[(ncbi_down_df.local_filename.str[-19:] == '_assembly_stats.txt') & # one row for every accession
                                (ncbi_down_df.organism_name.str.find(' genomosp.') == -1) & # skip not annotated species <1>
                                (ncbi_down_df.organism_name.str.find(' genosp.') == -1) & # skip not annotated species <2>
                                (ncbi_down_df.organism_name.str.find(' sp.') == -1) ] # skip not annotated species <3>

    print(ncbi_down_df)
    table = ncbi_down_df.to_csv('all_refseq_bacteria_unique.tsv', index=False, header=True, sep='\t')
