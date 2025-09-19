library("optparse")

parser <- OptionParser()
parser <- add_option(parser, c("-c", "--counts"), type="character", 
                    help="Path to cluster counts file")
parser <- add_option(parser, c("-f", "--filtered_counts"), type="character", 
                    help="Path to filtered cluster counts file")
parser <- add_option(parser, c("-t", "--taxonomy"), type="character", 
                    help="Path to cluster taxonomy file")
parser <- add_option(parser, c("-m", "--metadata"), type="character",
                    help="Path to metadata file")
parser <- add_option(parser, c("--sample_type_column"), type="character", default="lab_sample_type",
                    help="Column in metadata file that contains sample type information (default: lab_sample_type)")
parser <- add_option(parser, c("--sample_types"), type="character", default="sample",
                    help="Comma-separated list of sample types in metadata (default: sample)")
parser <- add_option(parser, c("--control_types"), type="character", default=paste0(c("buffer_blank","pcr_neg","extraction_neg","buffer_blank_art_spikes"), collapse=","),
                    help="Comma-separated list of control types in metadata (default: buffer_blank,pcr_neg,extraction_neg,buffer_blank_art_spikes)")
parser <- add_option(parser, c("--control_cutoff"), type="numeric", default=0.05,
                    help="Threshold for removing control clusters. Clusters occurring in more than control_cutoff of control samples will be removed. (default: 0.05)")
parser <- add_option(parser, c("--spikein_cutoff"), type="numeric", default=0.75,
                    help="Threshold for identifying spikein clusters. Clusters occurring in more than spikein_cutoff of spikein samples will be identified as spikeins. (default: 0.8)")
parser <- add_option(parser, c("--ignore_spikes"), action="store_true", default=FALSE,
                    help="Ignore spikein samples")
parser <- add_option(parser, c("--skip_control_cleaning"), action="store_true", default=FALSE,
                    help="Skip identification and removal of control clusters")                    
parser <- add_option(parser, c("--spikein_column"), type="character", default="biological_spikes",
                    help="Name of column in metadata file that determines whether a sample is a spikein (default: biological_spikes)")
parser <- add_option(parser, c("--counts_outfile"), type="character", default="cleaned_filtered_counts.tsv",
                    help="Path to cleaned filtered counts output file (default: cleaned_filtered_counts.tsv)")
parser <- add_option(parser, c("--taxonomy_outfile"), type="character", default="cleaned_filtered_cluster_taxonomy.tsv",
                    help="Path to cleaned filtered taxonomy output file (default: cleaned_filtered_cluster_taxonomy.tsv)")
parser <- add_option(parser, c("--control_outfile"), type="character", default=NULL,
                    help="Path to removed control clusters file. If specified, control clusters will be written to this file.")
parser <- add_option(parser, c("--spikein_outfile"), type="character", default=NULL,
                    help="Path to removed spikein clusters file. If specified, spikein clusters will be written to this file.")
parser <- add_option(parser, c("-r", "--remove_taxa"), type="character",
                    help="Commad-separated list of rank:taxa combinations to remove")

args <- parse_args(parser)
args$sample_types <- unlist(strsplit(args$sample_types,","))
args$control_types <- unlist(strsplit(args$control_types,","))
args$remove_taxa <- unlist(strsplit(args$remove_taxa,","))
if (any(is.null(args$counts), is.null(args$filtered_counts), is.null(args$taxonomy), is.null(args$metadata))) {
    stop("Please provide paths to counts (-c), filtered counts (-f), taxonomy (-t), and metadata (-m) files. See --help for more information.")    
}

library(data.table)
source("spikes_controls_fxns.R")

cat(paste0("Reading in metadata from ", args$metadata, "\n"))
meta <- read.delim(args$metadata, row.names=1)

# Get samples and controls
if (args$sample_type_column %in% colnames(meta)) {
    samples <- rownames(meta)[meta[,args$sample_type_column] %in% args$sample_types]
    controls <- rownames(meta)[meta[, args$sample_type_column] %in% args$control_types]
} else {
    stop(paste0("Column ", args$sample_type_column, " not found in metadata file."))
}

# Get spikein samples
if (args$spikein_column %in% colnames(meta)) {
    meta$spikein_sample <- as.logical(meta[, args$spikein_column])
    spikein_samples <- rownames(meta)[meta$spikein_sample==TRUE]
} else {
    spikein_samples <- NULL
}

cat("Samples: ", length(samples), "\n")
cat("Controls: ", length(controls), "\n")
cat("Spikein samples: ", length(spikein_samples), "\n")

cat(paste0("Reading in filtered counts from ", args$filtered_counts, "\n"))
filtered_counts <- fread(args$filtered_counts)

cat(paste0("Reading in raw cluster counts from ", args$counts, "\n"))
raw_counts <- fread(args$counts)
idx <- which(colSums(raw_counts[,2:ncol(raw_counts)])!=0)
idx <- c(1,idx+1)
raw_counts <- raw_counts[,..idx]

cat(paste0("Reading in taxonomy from ", args$taxonomy, "\n"))
taxonomy <- read.delim(args$taxonomy)
#taxonomy <- taxonomy[taxonomy$representative==1,]

if (length(spikein_samples)>0) {
    cat("Identifying spikeins\n")
    res <- remove_spikes(filtered_counts,spikein_samples, taxonomy, args$spikein_cutoff)
    if (args$ignore_spikes) {
        cat("Skipping removal of spikeins\n")
        cleaned_filtered_counts <- filtered_counts
    } else {
        cleaned_filtered_counts <- res$counts
    }
    spikein_remove_tax <- res$spike_tax
} else {
    cat("No spikein-samples found\n")
    cleaned_filtered_counts <- filtered_counts
    spikein_remove_tax <- list()
    spikein_remove_tax$cluster <- c()
}

if (!args$skip_control_cleaning) {
    cat("Identifying control clusters\n")
    res <- remove_control_clusters(cleaned_filtered_counts, raw_counts,taxonomy, samples,controls, spikein_remove_tax$cluster, args$ignore_spikes, args$control_cutoff)
    cleaned_filtered_counts <- res$counts
    control_remove_tax <- res$remove_tax
} else {
    cat("Skipping identification of control clusters\n")
    control_remove_tax <- NULL
}

cleaned_filtered_taxonomy <- taxonomy[taxonomy$cluster %in% cleaned_filtered_counts$cluster,]

if (!is.null(args$remove_taxa)) {
    n <- nrow(cleaned_filtered_taxonomy)
    for (pair in args$remove_taxa) {
        cat(paste0("Removing ", pair, "\n"))
        pair <- unlist(strsplit(pair,":"))
        rank <- pair[1]
        taxa <- pair[2]
        cleaned_filtered_taxonomy <- cleaned_filtered_taxonomy[cleaned_filtered_taxonomy[rank]!=taxa, ]
        cat(paste0("Removed ", n-nrow(cleaned_filtered_taxonomy), "\n"))
        n <- nrow(cleaned_filtered_taxonomy)
    }
}

cleaned_filtered_counts <- cleaned_filtered_counts[cluster%in%cleaned_filtered_taxonomy$cluster,]

cat(paste0("Writing cleaned filtered counts to ", args$counts_outfile, "\n"))
write.table(cleaned_filtered_counts,args$counts_outfile,sep="\t",quote=FALSE,row.names=FALSE)

cat(paste0("Writing cleaned filtered cluster taxonomy to ", args$taxonomy_outfile, "\n"))
write.table(cleaned_filtered_taxonomy,args$taxonomy_outfile,sep="\t",quote=FALSE,row.names=FALSE)

if ((!is.null(args$control_outfile)) & (!args$skip_control_cleaning)) {
    cat(paste0("Writing removed control clusters to ", args$control_outfile, "\n"))
    write.table(control_remove_tax,args$control_outfile,sep="\t",quote=FALSE,row.names=FALSE)
}

if (!is.null(args$spikein_outfile)) {
    cat(paste0("Writing removed spikein clusters to ", args$spikein_outfile, "\n"))
    write.table(spikein_remove_tax,args$spikein_outfile,sep="\t",quote=FALSE,row.names=FALSE)
}