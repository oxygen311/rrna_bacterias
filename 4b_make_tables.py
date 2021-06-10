import pandas as pd

import json

table = 'all_refseq_bacteria_unique.tsv'

skipped_16s, skipped_ori_ter = 0, 0


def get_species_subsp_strain(organism_name, infraspecific_name, isolate):
    organism_name = organism_name.replace('[', '').replace(']', '').replace('Candidatus ', '')

    subsp_splitted = organism_name.split(' subsp. ')
    subsp = subsp_splitted[1].split(' str. ')[0] if len(subsp_splitted) > 1 else ''

    organism_name_splitted = subsp_splitted[0].split(maxsplit=2)

    strain = organism_name_splitted[2] if len(organism_name_splitted) > 2 else ''
    if isinstance(infraspecific_name, str):
        infraspecific_name = infraspecific_name.replace('strain=', '')
        if strain in infraspecific_name: strain = infraspecific_name
        if not infraspecific_name in strain: strain += ('' if strain == '' else ' ') + infraspecific_name

    if isinstance(isolate, str):
        if not isolate in strain: strain += ('' if strain == '' else ' ') + isolate

    species = organism_name_splitted[0] + ' ' + organism_name_splitted[1]
    # print(f'<{species}>', f'<{subsp}>', f'<{strain}>', '\t\t', organism_name, '|', infraspecific_name, '|', isolate)

    return species, subsp, strain


def get_ori_ter(accession, asm):
    global skipped_ori_ter
    ori_ter_file = f'refseq/bacteria/{accession}/{accession}_{asm}_ori_ter.json'

    ans = []
    try:
        with open(ori_ter_file) as f:
            ans.append([accession, *json.load(f)])
    except FileNotFoundError:
        print('Skip (ori/ter):', accession)
        skipped_ori_ter += 1

    return ans


def get_16s(accession, asm):
    global skipped_16s
    ori_ter_file = f'refseq/bacteria/{accession}/{accession}_{asm}_16S.json'

    ans = []
    try:
        with open(ori_ter_file) as f:
            for j in json.load(f):
                ans.append([accession, *j])
    except FileNotFoundError:
        print('Skip (16s):', accession)
        skipped_16s += 1

    return ans


if __name__ == "__main__":
    ncbi_down_df = pd.read_csv(table, sep='\t')

    species_strain_2d = [[_1, *get_species_subsp_strain(_2, _3, _4)] for _1, _2, _3, _4 in
                         zip(ncbi_down_df.assembly_accession, ncbi_down_df.organism_name,
                             ncbi_down_df.infraspecific_name, ncbi_down_df.isolate)]
    ori_ter_2d = [result for _1, _2 in zip(ncbi_down_df.assembly_accession, ncbi_down_df.asm_name) for result in
                  get_ori_ter(_1, _2)]
    _16s_2d = [result for _1, _2 in zip(ncbi_down_df.assembly_accession, ncbi_down_df.asm_name) for result in
               get_16s(_1, _2)]

    species_strain_df = pd.DataFrame(data=species_strain_2d,
                                     columns=['assembly_accession', 'species', 'subsp', 'strain'])
    ori_ter_df = pd.DataFrame(data=ori_ter_2d,
                              columns=['assembly_accession', 'seq_id', 'origin', 'terminus', 'repl_1_mean',
                                       'repl_2_mean', 'c_skew_min', 'c_skew_max', 'length'])
    _16s_df = pd.DataFrame(data=_16s_2d,
                           columns=['assembly_accession', 'start', 'end', 'orientation'])

    species_strain_df.to_csv('species_subsp_strain.tsv', index=False, header=True, sep='\t')
    ori_ter_df.to_csv('ori_ter.tsv', index=False, header=True, sep='\t')
    _16s_df.to_csv('16S.tsv', index=False, header=True, sep='\t')

    print('Skipped 16s:', skipped_16s)
    print('Skipped ori/ter:', skipped_ori_ter)

    # table = ncbi_down_df.to_csv('ori_ter_table.tsv', index=False, header=True, sep='\t')
