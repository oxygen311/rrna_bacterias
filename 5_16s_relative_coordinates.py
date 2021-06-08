import pandas as pd


def relative_coordinates(start, end, origin, terminus, length):
    def get_relative_position(loc, s_terminus):
        return loc / s_terminus if loc < s_terminus else - (length - loc) / (length - s_terminus)

    s_start = (start - origin + length) % length
    s_end = (end - origin + length) % length
    s_terminus = (terminus - origin + length) % length

    rel_start = get_relative_position(s_start, s_terminus)
    rel_end = get_relative_position(s_end, s_terminus)

    return rel_start, rel_end


if __name__ == "__main__":
    df_16s = pd.read_csv('16S.tsv', sep='\t')
    df_ori_ter = pd.read_csv('ori_ter.tsv', sep='\t')
    df = pd.merge(df_16s, df_ori_ter, on='assembly_accession')

    print(df_16s)

    rel_coords_2d = [[_1, *relative_coordinates(_2, _3, _5, _6, _7), _4] for _1, _2, _3, _4, _5, _6, _7 in
                     zip(df.assembly_accession, df.start, df.end, df.orientation, df.origin, df.terminus, df.length)]

    rel_coords_df = pd.DataFrame(data=rel_coords_2d,
                                 columns=['assembly_accession', 'rel_start', 'rel_end', 'orientation'])

    rel_coords_df.to_csv('16S_rel.tsv', index=False, header=True, sep='\t')
