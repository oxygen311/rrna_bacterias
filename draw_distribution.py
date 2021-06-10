import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


bins_all = np.arange(-1, 1+1e-6, 0.025)
bins_grp = np.arange(-1, 1+1e-6, 0.0125)


def histogram_all():
    df_16s = pd.read_csv('16S_rel.tsv', sep='\t')
    df_16s['rel_center'] = (df_16s['rel_start'] + df_16s['rel_end']) / 2

    sns.histplot(data=df_16s, x='rel_center', bins=bins_all)

    plt.tight_layout()
    plt.xlim(xmin=-1, xmax=1)

    plt.savefig(f'charts/16s_relative_all.pdf')
    plt.show()


def histogram_grouped_species():
    df_16s = pd.read_csv('16S_rel.tsv', sep='\t')
    df_16s['rel_center'] = (df_16s['rel_start'] + df_16s['rel_end']).abs() / 2
    df_sp = pd.read_csv('species_subsp_strain.tsv', sep='\t')

    df = pd.merge(df_16s, df_sp, on='assembly_accession')

    for species, df_species in df.groupby('species'):
        strains = len(df_species.assembly_accession.unique())
        if strains < 50: continue
        print(species, len(df_species))

        plt.figure()
        sns.set_style('whitegrid')
        plt.title(f'{species}, {strains} strains, {round(len(df_species) / strains, 2)} mean 16s copies')
        sns.histplot(data=df_species, x='rel_center', bins=bins_grp)

        plt.tight_layout()
        plt.xlim(xmin=0, xmax=1)

        plt.savefig(f'charts/16s_relative_by_species/{"_".join(species.split())}_{strains}_strains.pdf')
        # plt.show()


def gc_skew():
    df = pd.read_csv('ori_ter.tsv', sep='\t')
    # df['skew'] = df['c_skew_max'] - df['c_skew_min']
    df['skew'] = (df['repl_1_mean'] - df['repl_2_mean']).abs()
    df['repr1_length'] = [(t - o + l) % l for o, t, l in zip(df.origin, df.terminus, df.length)]
    df['repr2_length'] = [(o - t + l) % l for o, t, l in zip(df.origin, df.terminus, df.length)]
    df['repr_length_diff'] = df['repr1_length'] / df['repr2_length']

    sns.histplot(data=df, x='skew')

    plt.tight_layout()
    plt.savefig(f'charts/mean_repl_skew.pdf')
    plt.show()


def dict_16s_copies(df_16s):
    return {asm: len(df_asm) for asm, df_asm in df_16s.groupby('assembly_accession')}


def annotate_16s_copies(df, _16s_copies):
    df['16s_copies'] = [_16s_copies.get(asm, -1) for asm in df.assembly_accession]


def scatter():
    df = pd.read_csv('ori_ter.tsv', sep='\t')

    # 16 copies
    _16s_copies = dict_16s_copies(pd.read_csv('16S_rel.tsv', sep='\t'))
    annotate_16s_copies(df, _16s_copies)

    # df['skew'] = df['c_skew_max'] - df['c_skew_min']
    df['skew'] = (df['repl_1_mean'] - df['repl_2_mean']).abs()
    df = df[df['16s_copies'] != -1]

    # df = df[(df.length < .8e7) & (df['skew'] < 0.3)]
    # sns.jointplot(data=df, x="16s_copies", y="skew", kind="hex")
    sns.boxplot(data=df, x="16s_copies", y="length")

    plt.tight_layout()
    plt.savefig(f'charts/boxplot_copies_length_filtered.pdf')
    plt.show()


def histogram_grouped_16s_copies():
    df_16s = pd.read_csv('16S_rel.tsv', sep='\t')
    annotate_16s_copies(df_16s, dict_16s_copies(df_16s))

    df_16s['rel_center'] = (df_16s['rel_start'] + df_16s['rel_end']).abs() / 2
    df_sp = pd.read_csv('species_subsp_strain.tsv', sep='\t')

    df = pd.merge(df_16s, df_sp, on='assembly_accession')

    for copies, df_copies in df.groupby('16s_copies'):
        strains = len(df_copies.assembly_accession.unique())
        if strains < 50: continue
        print(copies, len(df_copies))

        plt.title(f'{copies} copies, {strains} strains')
        sns.histplot(data=df_copies, x='rel_center', bins=bins_grp)

        plt.tight_layout()
        plt.xlim(xmin=0, xmax=1)

        plt.savefig(f'charts/16s_relative_by_16s_copies/{str(copies).zfill(2)}_copies_{strains}_strains_.pdf')
        plt.show()


def leading_lagging_strand():
    df_16s = pd.read_csv('16S_rel.tsv', sep='\t')
    df_16s['rel_center'] = (df_16s['rel_start'] + df_16s['rel_end']) / 2
    df_16s['rel_center_abs'] = df_16s['rel_center'].abs()

    df_16s['strand'] = ['leading' if (l > 0 and o == '+') or (l < 0 and o == '-') else 'lagging'
                      for l, o in zip(df_16s.rel_center, df_16s.orientation)]

    sns.histplot(data=df_16s, x='rel_center_abs', bins=bins_grp, hue='strand', stat='density',
                 element="step", common_norm=False)

    # sns.histplot(df_16s[df_16s.strand == 'leading'], x='rel_center_abs', bins=bins_grp, alpha=0.4, color='orange',
    #              label='leading strand', stat="density")
    # sns.histplot(df_16s[df_16s.strand != 'leading'], x='rel_center_abs', bins=bins_grp, alpha=0.3,
    #              label='lagging strand', stat="density")
    #
    plt.tight_layout()
    plt.xlim(xmin=0, xmax=1)

    plt.savefig(f'charts/leading_lagging_strand_density.pdf')
    plt.show()


if __name__ == "__main__":
    sns.set_style('whitegrid')
    leading_lagging_strand()