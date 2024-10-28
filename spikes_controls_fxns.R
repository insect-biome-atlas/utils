library(data.table)


# Help function to get prop_samples, mean, and max reads
#   index:  Index of rows to remove (TRUE) or keep (FALSE)
mean_max <- function(dt,index) {

    is_missing <- rowSums(dt[,2:ncol(dt)])==0

    remove_stats <- index & !is_missing
    remove_mean <- rowSums(dt[remove_stats,2:ncol(dt)]) / rowSums(dt[remove_stats,2:ncol(dt)]>0)
    remove_max <- apply(as.matrix(dt[remove_stats,2:ncol(dt)]),1,max)
    remove_prop <- rowMeans(dt[remove_stats,2:ncol(dt)]>0)

    is_missing_keep <- rowSums(dt[!index,2:ncol(dt)])==0
    keep_stats <- !index & !is_missing
    keep_mean <- rowSums(dt[keep_stats,2:ncol(dt)]) / rowSums(dt[keep_stats,2:ncol(dt)]>0)
    keep_max <- apply(as.matrix(dt[keep_stats,2:ncol(dt)]),1,max)
    keep_prop <- rowMeans(dt[keep_stats,2:ncol(dt)]>0)

    list(remove_mean=remove_mean,
         remove_max=remove_max,
         remove_prop=remove_prop,
         remove_missing=sum(is_missing & index),
         keep_mean=keep_mean,
         keep_max=keep_max,
         keep_prop=keep_prop,
         keep_missing=sum(is_missing & !index))
}


# Function to identify control clusters
#
#   counts:     cluster counts for all samples (including controls)
#   samples:    names of samples
#   controls:   names of controls
#   cutoff:     threshold for removing control clusters
identify_control_clusters <- function(counts, taxonomy, samples, controls, cutoff=0.05) {

    # Identify clusters that appear in controls above cutoff
    orig_tax <- taxonomy
    taxonomy <- taxonomy[taxonomy$representative==1,]
    # Extract sample and control reads
    idx <- which(colnames(counts) %in% samples)
    idx <- c(1,idx) # keep cluster name
    sample_counts <- counts[,..idx]

    idx <- which(colnames(counts) %in% controls)
    idx <- c(1,idx) # keep cluster name
    control_counts <- counts[,..idx]

    # identify clusters that have at least one read in a control
    include_rows <- rowSums(control_counts[, 2:ncol(control_counts)])>0
    # get counts of the control clusters in the controls
    control_counts <- control_counts[include_rows,]
    # get counts of the control clusters in the samples
    sample_counts <- sample_counts[include_rows,]

    # Output some diagnostics
    tot_clusters <- nrow(control_counts)
    # prop_controls is the proportion of control samples in which each cluster occurs
    prop_controls <- rowMeans(control_counts[,2:ncol(control_counts)]>0)

    cat("There are", tot_clusters, "clusters in the control samples\n")
    cat("Of these,", sum(prop_controls>cutoff), "occur in more than", cutoff, "of samples\n")
    # identify clusters that occur in more than cutoff proportion of control samples
    remove_clusters <- control_counts$cluster[prop_controls>cutoff]
    
    res <- mean_max(control_counts, prop_controls>cutoff)
    cat("Reads in controls of removed clusters:\n")
    remove_tax <- data.frame(list(cluster=remove_clusters,
                                  prop_controls=res$remove_prop,
                                  mean_reads=res$remove_mean,
                                  max_reads=res$remove_max,
                                  Genus=taxonomy$Genus[match(remove_clusters,taxonomy$cluster)],
                                  Species=taxonomy$Species[match(remove_clusters,taxonomy$cluster)],
                                  BOLD_bin=taxonomy$BOLD_bin[match(remove_clusters,taxonomy$cluster)]))
    #print(remove_tax)
    cat("Summary of reads in controls of kept clusters:\n")
    cat("prop_samples:\n")
    print(summary(res$keep_prop))
    cat("mean:\n")
    print(summary(res$keep_mean))
    cat("max:\n")
    print(summary(res$keep_max))

    prop_samples <- rowMeans(sample_counts[,2:ncol(sample_counts)]>0)
    res <- mean_max(sample_counts, prop_controls>cutoff)
    cat("Reads in samples of removed clusters:\n")
    remove_tax <- data.frame(list(cluster=remove_clusters,
                                  prop_samples=res$remove_prop,
                                  mean_reads=res$remove_mean,
                                  max_reads=res$remove_max,
                                  Genus=taxonomy$Genus[match(remove_clusters,taxonomy$cluster)],
                                  Species=taxonomy$Species[match(remove_clusters,taxonomy$cluster)],
                                  BOLD_bin=taxonomy$BOLD_bin[match(remove_clusters,taxonomy$cluster)]))
    print(remove_tax)
    cat("Summary of reads in samples of kept clusters:\n")
    cat("prop_samples:\n")
    print(summary(res$keep_prop))
    cat("mean:\n")
    print(summary(res$keep_mean))
    cat("max:\n")
    print(summary(res$keep_max))

    list(remove_clusters=remove_clusters, remove_tax=remove_tax)
}


# Function to remove control_clusters
#
#   filtered_counts:    filtered cluster counts for samples, col 1 should be cluster name
#   all_cluster_counts: cluster counts for all samples
#   samples:            names of samples
#   controls:           names of controls
#   cutoff:             threshold for removing control clusters
remove_control_clusters <- function(filtered_counts, all_cluster_counts, taxonomy, samples, controls, cutoff=0.05) {

    res <- identify_control_clusters(all_cluster_counts, taxonomy, samples, controls, cutoff)
    remove_clusters <- res$remove_clusters
    remove_tax <- res$remove_tax

    res <- list(
        counts=filtered_counts[!filtered_counts$cluster %in% remove_clusters,], 
        filtered_tax=taxonomy[!taxonomy$cluster %in% remove_clusters,], 
        remove_tax=remove_tax)
    res
}


# Function for identifying spikeins
identify_spikes <- function(counts, spikein_samples, taxonomy, cutoff=0.8) {

    # Get a rep asv taxonomy in case a complete cluster taxonomy is provided
    if ("representative" %in% colnames(taxonomy))
        taxonomy <- taxonomy[taxonomy$representative==1,]

    # Get counts for the samples containing spikeins
    # Account for counts being either data.frame or data.table
    idx <- which(colnames(counts) %in% spikein_samples)
    idx <- c(1,idx)
    if (class(counts)[1]=="data.table")
        counts <- counts[,..idx]
     else
        counts <- counts[,idx]

    # Identify spikeins
    prop_samples <- rowMeans(counts[,2:ncol(counts)]>0)
    spikein_candidates <- counts$cluster[prop_samples>cutoff]
    synthetic_spikeins <- spikein_candidates[is.na(match(spikein_candidates,taxonomy$cluster))]
    if (length(synthetic_spikeins>0)) {
        cat("Found", length(synthetic_spikeins), "synthetic spikein clusters: ")
        cat(synthetic_spikeins,sep=",")
        cat("\n")
    }
    bio_spikein_candidates <- spikein_candidates[!(spikein_candidates %in% synthetic_spikeins)]
    bio_spikeins <- bio_spikein_candidates[taxonomy$Class[match(bio_spikein_candidates,taxonomy$cluster)]=="Insecta"]
    cat("Found", length(bio_spikeins), "biological spikein clusters:\n")
    spike_tax <- taxonomy[taxonomy$cluster %in% bio_spikeins,c("cluster","Genus","Species","BOLD_bin")]
    spike_tax$prop_samples <- rowMeans(counts[match(spike_tax$cluster,counts$cluster),2:ncol(counts)]>0)
    print(spike_tax)

    # Return clusters
    res <- list(spikein_clusters=c(synthetic_spikeins, bio_spikeins), spike_tax=spike_tax)
    res
}


# Function for removing spikes
remove_spikes <- function(counts, spikein_samples, taxonomy, cutoff=0.8) {

    res <- identify_spikes(counts, spikein_samples, taxonomy, cutoff)
    spikein_clusters <- res$spikein_clusters
    spike_tax <- res$spike_tax

    res <- list(counts=counts[!(counts$cluster %in% spikein_clusters),], spike_tax=spike_tax)
    res
}

