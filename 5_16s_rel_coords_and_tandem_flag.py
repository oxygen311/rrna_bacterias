import pandas as pd
import numpy as np


def relative_coordinates(start, end, origin, terminus, length):
    def get_relative_position(loc, s_terminus):
        return loc / s_terminus if loc < s_terminus else - (length - loc) / (length - s_terminus)

    s_start = (start - origin + length) % length
    s_end = (end - origin + length) % length
    s_terminus = (terminus - origin + length) % length

    rel_start = get_relative_position(s_start, s_terminus)
    rel_end = get_relative_position(s_end, s_terminus)

    return rel_start, rel_end


def is_tandem_annotation(df, threshhold=1e4):
    aa_pos_to_tandem = {}
    df['center'] = (df.start + df.end) / 2

    for asm, df_asm in df.groupby('assembly_accession'):
        ps = df_asm.center.sort_values().to_numpy()
        length = df_asm.length.iloc[0]

        if len(ps) > 1:
            ds_to_prev = (ps - np.roll(ps, 1) + length) % length
            ds_to_next = (np.roll(ps, -1) - ps + length) % length

            tandem_flags = np.logical_or(ds_to_prev < threshhold, ds_to_next < threshhold)

            for p, tandem_flag in zip(ps, tandem_flags):
                aa_pos_to_tandem[(asm, p)] = tandem_flag

    df['has_tandem'] = [aa_pos_to_tandem.get((asm, p), False) for asm, p in zip(df.assembly_accession, df.center)]


if __name__ == "__main__":
    df_16s = pd.read_csv('16S.tsv', sep='\t')

    df_ori_ter = pd.read_csv('ori_ter.tsv', sep='\t')
    df = pd.merge(df_16s, df_ori_ter, on='assembly_accession')

    is_tandem_annotation(df)

    rel_coords_2d = [[_1, *relative_coordinates(_2, _3, _5, _6, _7), _4, _8] for _1, _2, _3, _4, _5, _6, _7, _8 in
                     zip(df.assembly_accession, df.start, df.end, df.orientation, df.origin, df.terminus, df.length,
                         df.has_tandem)]

    rel_coords_df = pd.DataFrame(data=rel_coords_2d,
                                 columns=['assembly_accession', 'rel_start', 'rel_end', 'orientation', 'has_tandem'])

    rel_coords_df.to_csv('16S_rel.tsv', index=False, header=True, sep='\t')
