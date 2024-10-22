# Insect Biome Atlas utils

R scripts useful for working with the IBA data.


## Installation


1. Clone the repository
```bash
git clone git@github.com:insect-biome-atlas/utils.git
cd utils
```

2. Install the required packages
```bash
conda env create
```

3. Activate the conda environment
```bash
conda activate iba-utils
```

## Usage

### Cleaning ASV data

The Rscript `clean_asv_data.R` can be used to clean the ASV data. The script can
remove ASV clusters present in control samples (_e.g._ negative controls, buffer
blanks etc.) as well as clusters that represent spike-ins. The script takes the
following arguments:

- `-c`, `--counts`: Path to cluster counts file

This should be a tab-separated file with ASV clusters as rows and samples
as columns. This represents the raw counts of ASV clusters prior to any
filtering.

Example:

| cluster      | sample1 | sample2 | sample3 |
|--------------|---------|---------|---------|
| ASV_cluster1 | 10      |       0 |       5 |
| ASV_cluster3 |  0      |       0 |      15 |
| ASV_cluster3 |  0      |      10 |       0 |

- `-f`, `--filtered_counts`: Path to filtered cluster counts file

This should be a tab-separated file with ASV clusters as rows and samples
as columns. This represents filtered counts of ASV clusters after removing
noise such as NUMTs and low abundance clusters.

Example:

| cluster      | sample1 | sample2 | sample3 |
|--------------|---------|---------|---------|
| ASV_cluster1 | 10      |       0 |       5 |
| ASV_cluster3 |  0      |       0 |      15 |
| ASV_cluster3 |  0      |      10 |       0 |

- `-t`, `--taxonomy`: Path to cluster taxonomy file

This should be a tab-separated file with ASV ids in the first column. The
file **must** contain a column `cluster` with ASV cluster designations as well
as a column `representative` with values either `1` or `0` indicating whether
the ASV is a representative of the cluster.

Example:

| ASV | cluster      | representative | Kingdom  | Phylum     | Class   | ... |
|-----|--------------|----------------|----------|------------|---------|-----|
| ASV1| ASV_cluster1 |              1 | Animalia | Arthropoda | Insecta | ... |
| ASV2| ASV_cluster1 |              0 | Animalia | Arthropoda | Insecta | ... |
| ASV3| ASV_cluster2 |              1 | Animalia | Arthropoda | Insecta | ... |
| ASV4| ASV_cluster3 |              1 | Animalia | Arthropoda | Insecta | ... |

- `-m`, `--metadata`: Path to metadata file

This should be a tab-separated file with sample ids in the first column. The
file **must** contain a column that designates the type of sample allowing the
script to discriminate between true samples and controls (see options
`--sample_type_column`, `--sample_types` and `--control_types` below). If
spike-in clusters are to be removed, there must also be a column named
`spikein_sample` that contains `1` or `True` for samples to which spike-ins were
added.

Example:

| sample | lab_sample_type | spikein_sample |
|--------|-----------------|----------------|
| sample1| sample          | 0              |
| sample2| buffer_blank    | 0              |
| sample3| pcr_neg         | 0              |
| sample4| sample          | 1              |


- `--sample_type_column`: Column in metadata file that contains sample type (default: `lab_sample_type`)

Name of column that contains the sample type. This is used to discriminate between
true samples and controls. The default value is `lab_sample_type`.

- `--sample_types`: Comma-separated list of sample types (default: `sample`)

Comma-separated list of sample types. These are the sample types that are considered
true samples and not controls. The default value is `sample`.

- `--control_types`: Comma-separated list of control types (default: `buffer_blank`, `pcr_neg`, `extraction_neg`, `buffer_blank_art_spikes`)

Comma-separated list of control types. These are the sample types that are considered
controls and which are used to remove ASV clusters based on occurrence in these
samples. The default value is `buffer_blank`, `pcr_neg`, `extraction_neg`, `buffer_blank_art_spikes`.

- `--control_cutoff`: Threshold for removing control clusters (default: 0.05)

Clusters occurring in more than `control_cutoff` of control samples will be
removed. The default value is `0.05`.

- `--spikein_cutoff`: Threshold for identifying spikein clusters (default: 0.8)

Clusters occurring in more than `spikein_cutoff` of spikein samples will be
identified as spikeins. The default value is `0.8`.

- `--counts_outfile`: Path to output file with counts of cleaned clusters (default: `cleaned_filtered_counts.tsv`)

This will be a tab-separated file with ASV clusters as rows and samples as
columns. This file will contain the filtered counts of ASV clusters after
removing noise, controls and spike-ins.

- `--taxonomy_outfile`: Path to output file with taxonomy of cleaned clusters (default: `cleaned_cluster_taxonomy.tsv`)

This will be a tab-separated file with ASV ids in the first column. The file
will contain the taxonomy of ASV clusters after removing noise, controls and
spike-ins.

- `--control_outfile`: Path to output file with clusters identified in control samples

This will be a tab-separated file with ASV clusters identified in control
samples in the first column and with additional columns containing read statistics and taxonomic assignments.

- `--spikein_outfile`: Path to output file with clusters identified as spike-ins

This will be a tab-separated file with ASV clusters identified as spike-ins in
the first column and with additional columns containing read statistics and taxonomic assignments.

Example of usage:
    
```bash
Rscript clean_asv_data.R -c cluster_counts.tsv \
    -f noise_filtered_cluster_counts.tsv \
    -t cluster_taxonomy.tsv \
    -m metadata.tsv \
    --counts_outfile cleaned_filtered_counts.tsv \
    --taxonomy_outfile cleaned_cluster_taxonomy.tsv
```