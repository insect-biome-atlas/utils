# Cleaning the IBA data

The commands below remove OTUs present in 5% of control samples, identifies
biological spike-ins by looking for OTUs present in at least 75% of samples (but
does not remove them) and removes the genus Zoarces from the data.

## Sweden

```bash
Rscript --vanilla clean_asv_data.R \
    -c data/cluster_counts_SE.tsv \
    -f data/noise_filtered_cluster_counts_SE.tsv \
    -t data/cluster_taxonomy_SE.tsv \
    -m data/CO1_sequencing_metadata_SE.tsv \
    --spikein_column biological_spikes --ignore_spikes \
    --counts_outfile results/cleaned_noise_filtered_cluster_counts_SE.tsv \
    --taxonomy_outfile results/cleaned_noise_filtered_cluster_taxonomy_SE.tsv \
     --control_outfile results/removed_control_tax_SE.tsv \
     --spikein_outfile results/spikeins_tax_SE.tsv \
     --remove_taxa Genus:Zoarces
```

## Madagascar

```bash
Rscript --vanilla clean_asv_data.R \
    -c data/cluster_counts_MG.tsv \
    -f data/noise_filtered_cluster_counts_MG.tsv \
    -t data/cluster_taxonomy_MG.tsv \
    -m data/CO1_sequencing_metadata_MG.tsv \
    --spikein_column biological_spikes --ignore_spikes \
    --counts_outfile results/cleaned_noise_filtered_cluster_counts_MG.tsv \
    --taxonomy_outfile results/cleaned_noise_filtered_cluster_taxonomy_MG.tsv \
     --control_outfile results/removed_control_tax_MG.tsv \
     --spikein_outfile results/spikeins_tax_MG.tsv \
     --remove_taxa Genus:Zoarces
```
