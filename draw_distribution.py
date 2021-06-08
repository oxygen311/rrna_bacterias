import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


bins_all = np.arange(-1, 1+1e-6, 0.025)
bins_grp = np.arange(-1, 1+1e-6, 0.05)

def histogram_all():
    df_16s = pd.read_csv('16S_rel.tsv', sep='\t')
    df_16s['rel_center'] = (df_16s['rel_start'] + df_16s['rel_end']) / 2

    sns.set_style('whitegrid')
    sns.histplot(data=df_16s, x='rel_center', bins=bins_all)

    plt.tight_layout()
    plt.xlim(xmin=-1, xmax=1)

    plt.savefig(f'charts/16s_relative_all.pdf')
    plt.show()

def histogram_grouped():
    df_16s = pd.read_csv('16S_rel.tsv', sep='\t')
    df_16s['rel_center'] = (df_16s['rel_start'] + df_16s['rel_end']) / 2
    df_sp = pd.read_csv('species_subsp_strain.tsv', sep='\t')

    df = pd.merge(df_16s, df_sp, on='assembly_accession')

    for species, df_species in df.groupby('species'):
        if len(df_species) < 50: continue
        print(species, len(df_species))

        sns.set_style('whitegrid')
        plt.title(f'{species}, {len(df_species)} strains')
        sns.histplot(data=df_species, x='rel_center', bins=bins_grp)

        plt.tight_layout()
        plt.xlim(xmin=-1, xmax=1)

        plt.savefig(f'charts/16s_relative_by_species/{"_".join(species.split())}_{len(df_species)}.pdf')
        plt.show()

if __name__ == "__main__":
    histogram_grouped()