# Code relies on data.table
library(data.table)

# Function for getting the IBA data for CO1.
#
# Parameters
#  -     data_path: path to the processed data
#  - metadata_path: path to the metadata
#  -       country: the country for which data is requested
#  -       dataset: the dataset requested; you can request more than one, but reads are not merged across datasets
#  -     calibrate: whether to calibrate the read counts using biological spike-in data (only possible for lysates and homogenates)
#  -                NB! when calibrate is requested, spike-ins will be removed
#  - remove_spikes: whether to remove reads from biological spike-ins (only relevant for lysates and homogenates)

get_iba_co1_data <- function(data_path,
                             metadata_path=data_path,
                             country=c("mg","se"),
                             dataset=c("homogenate|lysate|ethanol|soil|litter"),
                             calibrate=FALSE,
                             remove_spikes=TRUE) {
    
    # Refuse to lump datasets across countries
    if (country!="mg" && country !="se") {
        cat ("ERROR: 'country' must be one of 'mg' or 'se'\n")
        return (NA)
    }

    # Get data files
    cat("reading data files, this may take a while...\n")
    counts <- fread(paste0(data_path,"non_numts_clusters_counts.cleaned.tsv"),header=TRUE,sep="\t")
    taxonomy <- fread(paste0(data_path,"cluster_consensus_taxonomy.tsv"))

    # Initialize index into columns we want to keep
    index <- numeric()

    if (country=="se") {

        # Get metadata file and remove non-samples and sequencing failures
        meta <- fread(paste0(metadata_path,"CO1_sequencing_metadata_SE.tsv"))
        meta <- meta[meta$lab_sample_type=="sample" & meta$sequencing_status=="sequencing successful",]

        # Get and possibly adjust malaise trap data. We use the fact here
        # that the spike-ins are the same in both cases.
        if (grepl("lysate", dataset) || grepl("homogenate",dataset)) {
            malaise_samples <- character()
            if (grepl("lysate", dataset))
                malaise_samples <- with(meta, sampleID_NGI[dataset=="CO1_lysate_2019_SE"])
            if (grepl("homogenate", dataset)) {
                homogenates <- with(meta, sampleID_NGI[dataset=="CO1_homogenate_2019_SE"])
                malaise_samples <- c(malaise_samples, homogenates)
            }
            counts <- handle_spikes(counts, malaise_samples, taxonomy, calibrate, remove_spikes)
            index <- c(index, match(malaise_samples,colnames(counts)))
        }

        # Get ethanol data if requested
        if (grepl("ethanol", dataset)) {

            ethanol_samples <- with(meta, sampleID_NGI[dataset=="CO1_ethanol_2019_SE"])
            index <- c(index, match(ethanol_samples,colnames(counts)))
        }
 
        # Get soil and/or litter data if requested
        if (grepl("soil", dataset)) {

            soil_samples <- with(meta, sampleID_NGI[grepl("_S",sampleID_FIELD)])
            index <- c(index, match(soil_samples,colnames(counts)))
        }
        if (grepl("litter", dataset)) {

            litter_samples <- with(meta, sampleID_NGI[grepl("_L",sampleID_FIELD) & lab_sample_type=="sample"])
            index <- c(index, match(litter_samples,colnames(counts)))
        }
    }
    if (country=="mg") {

        # Get metadata file and remove non-samples and sequencing failures
        # Note misspelling of sequencing status value in MG metadata file
        # We are ignoring the single sample with very few reads
        meta <- fread(paste0(metadata_path,"CO1_sequencing_metadata_MG.tsv"))
        meta <- meta[meta$lab_sample_type=="sample" & meta$sequencing_status=="sequencing successfull",]

        # Adjust lysate data if requested
        if (grepl("lysate", dataset)) {

            lysates <- with(meta, sampleID_NGI[dataset=="CO1_lysate_2019_MG"])
            counts <- handle_spikes(counts, lysates, taxonomy, calibrate, remove_spikes)
            index <- c(index, match(lysates,colnames(counts)))
        }

        # Get litter data if requested
        if (grepl("litter", dataset)) {

            litter_samples <- with(meta, sampleID_NGI[dataset=="CO1_litter_2019_MG"])
            index <- c(index, match(litter_samples,colnames(counts)))
        }
    }
 
    # Extract indices, do not forget to keep the cluster column. Also, safeguard against
    # samples that have no counts data despite the metadata claiming the contrary
    index <- c(1, index)
    index <- index[!is.na(index)]
    counts <- counts[,..index]

    # Remove clusters that are not encountered in this dataset
    tot <- rowSums(counts[,2:ncol(counts)])
    counts <- counts[tot!=0,]

    # Add in taxonomy data and convert to data.frame
    data.frame(merge(counts, taxonomy, by="cluster"))
}


# Help function for calibration and removal of spike-ins
# We identify spike-in clusters as those that occur in
# more than 80% of samples and belong to Insecta
handle_spikes <- function(counts, samples, taxonomy, calibrate, remove_spikes) {

    # get column indices in the counts data table
    # safe_guard against samples that are missing in the counts table
    idx <- match(samples,colnames(counts))
    idx <- idx[!is.na(idx)]

    # identify spikeins
    spikein_candidates <- counts$cluster[rowMeans(counts[,..idx]>0)>0.8]
    spikein_clusters <- spikein_candidates[taxonomy$Class[match(spikein_candidates,taxonomy$cluster)]=="Insecta"]
    cat("found", length(spikein_clusters), "spikein clusters:\n")
    print(spikein_clusters)

    # calibrate
    # This is slow on data.table objects so it makes sense to convert to a matrix
    # and compute on the matrix; note that index differs by 1 because cluster column
    # is not retained in the matrix.
    # Also note that there are occasional samples without spike-ins; we simply do
    # not correct these read counts (what else can we do?). Presumably, spike-ins
    # were not added to these samples by mistake.
    if (calibrate && length(spikein_clusters) > 0) {
        cat("calibrating, this may take a while...\n")
        spike_counts <- colSums(counts[counts$cluster %in% spikein_clusters,..idx])
        correction <- spike_counts / mean(spike_counts[spike_counts!=0])
        m <- as.matrix(counts[,2:ncol(counts)])
        for (i in 1:length(idx)) {
            if (spike_counts[i]!=0)
                m[,idx[i]-1] <- ceiling(m[,idx[i]-1] / correction[i])
        }

        counts <- cbind(counts[,1],as.data.table(m))
    }

    # Remove the spike-ins
    if (calibrate || remove_spikes) {
        counts <- counts[!(counts$cluster %in% spikein_clusters),]
    }

    counts
} 

