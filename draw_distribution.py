import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

from collections import Counter


bins_all = np.arange(-1, 1+1e-6, 0.0125)
bins_grp = np.arange(-1, 1+1e-6, 0.025)


def histogram_all():
    df_16s = pd.read_csv('16S_rel.tsv', sep='\t')
    df_16s['rel_center'] = (df_16s['rel_start'] + df_16s['rel_end']).abs() / 2

    sns.histplot(data=df_16s, x='rel_center', bins=bins_all)

    plt.tight_layout()
    plt.xlim(xmin=0, xmax=1)

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


def dict_16s_copies(df_16s):
    return {asm: len(df_asm) for asm, df_asm in df_16s.groupby('assembly_accession')}


def annotate_16s_copies(df, _16s_copies):
    df['16s_copies'] = [_16s_copies.get(asm, -1) for asm in df.assembly_accession]


def boxplot():
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

    plt.title(Counter(df_16s.strand))

    plt.tight_layout()
    plt.xlim(xmin=0, xmax=1)

    plt.savefig(f'charts/leading_lagging_strand_density.pdf')
    plt.show()


def rank_copies_histogram():
    def rank_copies_from_origin():
        def make_rank(aas):
            curr_aa, rank, ans = '', 0, []
            for aa in aas:
                if aa == curr_aa:
                    rank += 1
                else:
                    curr_aa = aa
                    rank = 1
                ans.append(rank)
            return ans

        df_16s = pd.read_csv('16S_rel.tsv', sep='\t')
        annotate_16s_copies(df_16s, dict_16s_copies(df_16s))
        df_16s['rel_center_abs'] = (df_16s['rel_start'] + df_16s['rel_end']).abs() / 2

        df_16s = df_16s.sort_values(['assembly_accession', 'rel_center_abs'])
        df_16s['rankk'] = make_rank(df_16s.assembly_accession)

        return df_16s

    df_16s = rank_copies_from_origin()
    for copies, df_copies in df_16s.groupby('16s_copies'):
        strains = len(df_copies.assembly_accession.unique())
        if strains < 50: continue
        print(copies, len(df_copies))

        plt.title(f'{copies} copies, {strains} strains')

        # sns.histplot(data=df_copies, x='rel_center_abs', bins=bins_grp, hue='rankk', stat='density',
        #                            element="step", common_norm=False, palette='Set1')
        # plt.xlim(xmin=0, xmax=1)

        sns.violinplot(data=df_copies, y='rel_center_abs', x='rankk', common_norm=False, palette='Set1')

        plt.tight_layout()

        plt.savefig(f'charts/rank_copies_by_16s_copies/{str(copies).zfill(2)}_copies_{strains}_strains_vionil.pdf')
        plt.show()


def min_dist_btw_copies():
    def min_dist():
        aa_to_min_dist = {}

        for asm, df_asm in df_16s.groupby('assembly_accession'):
            ps = df_asm.rel_center.sort_values().to_numpy()

            if len(ps) > 1:
                ds = (ps - np.roll(ps, 1) + 1) % 1
                aa_to_min_dist[asm] = ds.min()

        return aa_to_min_dist

    df_16s = pd.read_csv('16S_rel.tsv', sep='\t')
    df_16s['rel_center'] = (df_16s['rel_start'] + df_16s['rel_end']) / 2

    aa_to_min_dist = min_dist()
    df = pd.read_csv('ori_ter.tsv', sep='\t')
    annotate_16s_copies(df, dict_16s_copies(df_16s))

    df['min_distance_btw_copies'] = [aa_to_min_dist.get(asm, -1) for asm in df.assembly_accession]
    df = df[df.min_distance_btw_copies != -1]

    for copies, df_copies in df.groupby('16s_copies'):
        strains = len(df_copies.assembly_accession.unique())
        if strains < 50: continue
        print(copies, len(df_copies))

        plt.title(f'{copies} copies, {strains} strains')

        sns.histplot(data=df_copies, x='min_distance_btw_copies', bins=bins_grp)
        plt.xlim(xmin=0, xmax=1)

        plt.tight_layout()

        plt.savefig(f'charts/minimal_distance_btw_copies_by_16s_copies/'
                    f'{str(copies).zfill(2)}_copies_{strains}_strains_hist.pdf')
        plt.show()


def has_tandem_distribution():
    df_16s = pd.read_csv('16S_rel.tsv', sep='\t')
    df_16s['rel_center'] = (df_16s['rel_start'] + df_16s['rel_end']).abs() / 2

    print(df_16s.columns)
    print(df_16s)

    # sns.histplot(data=df_16s, x='rel_center', bins=bins_grp, hue='has_tandem', stat='density',
    #              element="step", common_norm=False)

    sns.histplot(data=df_16s, x='rel_center', bins=bins_grp, hue='has_tandem', element="step")

    plt.title(Counter(df_16s.has_tandem))

    plt.tight_layout()
    plt.xlim(xmin=0, xmax=1)

    plt.savefig(f'charts/has_tandem_hist.pdf')
    plt.show()


def symmetry_2_copies(rand_size=1000000):
    df_16s = pd.read_csv('16S_rel.tsv', sep='\t')
    annotate_16s_copies(df_16s, dict_16s_copies(df_16s))

    df_16s['rel_center'] = (df_16s['rel_start'] + df_16s['rel_end']) / 2

    df_16s = df_16s[df_16s['16s_copies'] == 2]

    ds = []

    for asm, df_asm in df_16s.groupby('assembly_accession'):
        assert len(df_asm) == 2

        r1, r2 = df_asm.rel_center
        if (r1 > 0 and r2 > 0) or (r1 < 0 and r2 < 0): continue

        r = abs(abs(r1) - abs(r2))
        ds.append(['real', r])

    rnd1 = np.random.rand(rand_size) * 2 - 1
    rnd2 = np.random.rand(rand_size) * 2 - 1
    for r1, r2 in zip(rnd1, rnd2):
        if (r1 > 0 and r2 > 0) or (r1 < 0 and r2 < 0): continue
        r = abs(abs(r1) - abs(r2))
        ds.append(['random', r])

    rel_coords_df = pd.DataFrame(data=ds,
                                 columns=['type', 'symmetry'])
    sns.histplot(data=rel_coords_df, x='symmetry', hue='type', stat='density',
                 element="step", common_norm=False)

    plt.tight_layout()
    plt.xlim(xmin=0, xmax=1)

    plt.savefig(f'charts/symmetry_2_copies_density_v2.pdf')
    plt.show()


def copies_repl1_repl2():
    df_16s = pd.read_csv('16S_rel.tsv', sep='\t')
    df_16s['rel_center'] = (df_16s['rel_start'] + df_16s['rel_end']) / 2

    repr_12_2d = []
    for asm, df_asm in df_16s.groupby('assembly_accession'):
        repl_1 = (df_asm.rel_center >= 0).sum()
        repl_2 = len(df_asm) - repl_1

        if repl_2 > repl_1: repl_1, repl_2 = repl_2, repl_1
        repr_12_2d.append([repl_1, repl_2])

    repr_12_df = pd.DataFrame(data=repr_12_2d, columns=['repl1', 'repl2'])

    repr_count_df = repr_12_df.groupby(['repl1', 'repl2']).size().reset_index(name='counts')
    # repr_count_df = repr_count_df[repr_count_df.counts > 50]
    repr_count_df = repr_count_df.pivot("repl1", "repl2", "counts")

    sns.heatmap(data=repr_count_df, cmap='rocket_r')

    plt.tight_layout()
    plt.savefig('charts/repl1_repl2_count_heatmap.pdf')
    plt.show()


if __name__ == "__main__":
    sns.set_style('whitegrid')
    symmetry_2_copies()