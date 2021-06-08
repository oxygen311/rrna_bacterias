# 1) Download data

All bacterial complete assemblies were downloaded from RefSeq database (on 4 June 2021) was downloaded using `ncbi-genome-download` tool:

```bash
ncbi-genome-download bacteria --assembly-levels complete,chromosome --section refseq --formats fasta,rna-fasta,assembly-report,assembly-stats --parallel 32 --metadata-table all_refseq_bacteria.tsv --human-readable
```

Slurm script and log are named `ncbi-genome-download` with extension `.sh`, `.out` and `.err`. 

# 2) Origin and terminus annotation with

Function `gc_skew` from iRep package was used for origin and terminus annotation.
Script `2_ori_ter_annotation_single.py` makes annotation for single genome using this library.
Full slurm script for making every genome annotation separate subtask:

```bash
#SBATCH --array=1-25352
DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" folders.txt)
python3 ori_ter_annotation_single.py -f refseq/bacteria/${DIR}
```

Where `folders.txt` is a list of all folders containing genomes.

# 3) Readable info about 16S locations and orientations

Extraction info about ribosomal operons from feature table.
Script `3_16s_annotation_single.py` extracts 16S info including locations of start and end and orientation based on **16S ribosomal RNA** in name.
Full slurm script for making every genome annotation separate subtask:

```bash
#SBATCH --array=1-25352
DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" folders.txt)
python3 3_16s_annotation_single.py -f refseq/bacteria/${DIR}
```

# 4) Overall table with essentials info

## 4a) Filtering table

Script `4a_filter_table.py` filters `all_refseq_bacteria.tsv` table for keeping only one row for every strain.
Also, this script drops strains with not annotated species (`sp.` in `organism_name`).

## 4b) Generating tables

Script `4b_make_tables.py` generates tables with necessary information.
* Table `species_subsp_strain.tsv` with info about species (genus + species), subspecies if exists and strain names;
* Table `ori_ter.tsv` with info about origin and terminus locations;
* Table `16s.tsv` with info about presented 16S rRNAs, every row represents one rRNA.

All tables are equipped with `assembly_accession` column for further use.


# 5) Relative coordinates

Script `5_16s_relative_coordinates.py` generates table `16S_rel.tsv` with relative coordinates of 16S rRNA instead of absolute in a following way:

![Relative coordinates](figs/rel_coords.svg)
