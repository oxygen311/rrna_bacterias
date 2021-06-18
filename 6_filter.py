import pandas as pd


before_folder = 'before_filtering/'
after_folder = 'after_filtering/'


def consteuct_allowed():
    df = pd.read_csv(f'{before_folder}/ori_ter.tsv', sep='\t')

    df['skeww'] = (df['repl_1_mean'] - df['repl_2_mean']).abs()
    df['repr1_length'] = [(t - o + l) % l for o, t, l in zip(df.origin, df.terminus, df.length)]
    df['repr2_length'] = [(o - t + l) % l for o, t, l in zip(df.origin, df.terminus, df.length)]
    df['repr_length_diff'] = (df['repr1_length'] - df['repr2_length']).abs() \
                             / ((df['repr1_length'] + df['repr2_length']) / 2)
    df = df[(df.skeww > 0.05) & (df.repr_length_diff < 0.15) & (df.length > 800000)]

    df_16s = pd.read_csv(f'{before_folder}/16S.tsv', sep='\t')

    set1 = set(df.assembly_accession)
    set2 = set(df_16s.assembly_accession)

    return set1 & set2


def filter_file(file, allowed_genomes):
    df = pd.read_csv(before_folder + file, sep='\t')
    df = df[df.assembly_accession.isin(allowed_genomes)]
    df.to_csv(after_folder + file, index=False, header=True, sep='\t')


allowed_genomes = consteuct_allowed()
filter_file('16S.tsv', allowed_genomes)
filter_file('16S_rel.tsv', allowed_genomes)
filter_file('all_refseq_bacteria.tsv', allowed_genomes)
filter_file('all_refseq_bacteria_unique.tsv', allowed_genomes)
filter_file('ori_ter.tsv', allowed_genomes)
filter_file('species_subsp_strain.tsv', allowed_genomes)
