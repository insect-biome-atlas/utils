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
                    help="Column in metadata file that contains sample type")
parser <- add_option(parser, c("--sample_types"), type="character", default="sample",
                    help="Comma-separated list of sample types")
parser <- add_option(parser, c("--control_types"), type="character", default=c("buffer_blank","pcr_neg","extraction_neg","buffer_blank_art_spikes"),
                    help="Comma-separated list of control types")
parser <- add_option(parser, c("-o", "--outfile"), type="character", default="cleaned_filtered_counts.tsv",
                    help="Path to output file")

args <- parse_args(parser)
args$sample_types <- unlist(strsplit(args$sample_types,","))
args$control_types <- unlist(strsplit(args$control_types,","))
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
if ("spikein_sample" %in% colnames(meta)) {
    meta$spikein_sample <- as.logical(meta$spikein_sample)
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
taxonomy <- taxonomy[taxonomy$representative==1,]

cleaned_filtered_counts <- remove_control_clusters(filtered_counts,raw_counts,taxonomy,samples,controls)

if (length(spikein_samples)>0) {
    cleaned_filtered_counts <- remove_spikes(cleaned_filtered_counts,spikein_samples, taxonomy)
}

cat(paste0("Writing cleaned filtered counts to ", args$outfile, "\n"))
write.table(cleaned_filtered_counts,args$outfile,sep="\t",quote=FALSE,row.names=FALSE)