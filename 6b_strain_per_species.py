import pandas as pd


before_folder = 'after_filtering/'
after_folder = 'strain_per_species/'


def consteuct_allowed():
    df1 = pd.read_csv(f'{before_folder}species_subsp_strain.tsv', sep='\t')
    df2 = pd.read_csv(f'{before_folder}all_refseq_bacteria_unique.tsv', sep='\t')

    df = pd.merge(df1, df2, on='assembly_accession')

    allowed_genomes = set()

    for species, df_species in df.groupby('species'):
        if len(df_species) < 10: continue

        df_species_reference = df_species[df_species.refseq_category == 'reference genome']
        if len(df_species_reference) > 0:
            allowed_genomes.add(df_species_reference.assembly_accession.values[0])
            continue

        df_species_representative = df_species[df_species.refseq_category == 'representative genome']
        if len(df_species_representative) > 0:
            allowed_genomes.add(df_species_representative.assembly_accession.values[0])
            continue

        allowed_genomes.add(df_species.assembly_accession.values[0])

    return allowed_genomes


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
