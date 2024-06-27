# format community data -----------------------------------------------------------------------
library(data.table)

# load community data -----------------------------------------------------------------------------------

# Requires downloading of data files from figshare
IBA_taxa         <- fread("consensus.SE.taxonomy.tsv")  

# Clusters by sample
IBA_clusters     <- fread("CO1_lysate_2019_SE.cleaned.clusters.counts.tsv") 

# Metadata
IBA_seq_meta     <- fread("CO1_sequencing_metadata_SE.tsv")  
IBA_malaise_meta <- fread("samples_metadata_malaise_SE_2019.tsv") 
IBA_trap_meta    <- fread("traps_metadata_SE_2019.tsv") 


# Spike in IDs 
IBA_spike_ins  <- fread("malaise_bio_spikeins_SE_2019.tsv")
spike_in_clust <- unique(IBA_spike_ins$cluster) # Unique ID of spike-in clusters

# tidy data ----------------------------------------------------------------------------------

# Assemble sequence data
taxa_counts_DT <-  # Merge the taxonomy with the cluster data using merge.data.table()
                   merge(IBA_taxa, IBA_clusters, by = "cluster", all = TRUE) |>
                   # Create a long table with cluster info by each sample
                   melt(id.vars = 1:9 , variable.name = "sampleID_NGI" , value.name = "read_count") |> 
                   # Filter out clusters with no reads in a sample ( "_" is the base pipe placeholder)
                   _[read_count > 0,] 

# Assemble the metadata
full_meta_DT <- merge(IBA_seq_meta , IBA_malaise_meta , all=TRUE) |> 
                merge(y = _ , x = IBA_trap_meta , all=TRUE) |> 
                # Decide on which relevant columns to keep 
                _[,c("trapID","duration_min", "sampleID_NGI" ,"lab_sample_type")] 

# Combine metadata and OTU data
OTU_DT <- merge(taxa_counts_DT , full_meta_DT, by = "sampleID_NGI") 


# Additional filters
OTU_filtered_DT <- # Remove non-samples
                    OTU_DT[lab_sample_type == "sample",] 
                    # Filter out spike-ins 
                    OTU_DT[!cluster %in% spike_in_clust,]


