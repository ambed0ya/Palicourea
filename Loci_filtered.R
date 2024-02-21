##### Script to filter the loci based on their percentage of length recovered
##### and on the percentage of samples in which they were recovered (FROM LEO-PAUL DAGALLIER)


# Load the necessary packages ---------------------------------------------
library(ggpubr)
library(tidyverse)
library(reshape2)
library(wesanderson)

# Set the paths and load the data -----------------------------------------
# PATH SETTING AND DATA LOADING SHOULD BE DONE PRIOR TO CALL THIS SCRIPT
# THIS SHOULD BE DONE FOLLOWING THE COMMENTED EXAMPLE BELOW
path_to_data <- "/PATH_TO_FOLDER_WITH_HYBPIPER_OUTPUT/"
filename = paste0(path_to_data, "/seq_lengths.tsv")
seq_lengths_raw <- read.table(filename, header = T, row.names = 1, sep = "\t", check.names = F)
seqpath_to_out <- "/PATH_TO_FOLDER_WITH_HYBPIPER_OUTPUT/"
limit_perc_length_wanted = c(0.1, 0.20, 0.5, 0.75)
limit_perc_nb_wanted = c(0.1, 0.20, 0.5, 0.75)

# Prepare the data frames -------------------------------------------------
# FIRST: sort the locus names in the raw data frame
seq_lengths_raw <- seq_lengths_raw[,order(names(genes_sequences_lengths_raw))]
# get the length of reference sequence for each locus
reference_length = as.numeric(seq_lengths_raw[which(row.names(genes_sequences_lengths_raw)== "MeanLength"),])
names(reference_length) <- names(seq_lengths_raw)

# Keep the sequences lengths only for the genes (that is: remove MeaLength from the dataframe)
seq_lengths = seq_lengths_raw[which(row.names(seq_lengths_raw) != "MeanLength"),]

# Table of the number of samples per locus --------------------------------
sequence_PA = ifelse(seq_lengths > 0, 1, 0) # draw a presence/absence table
tmp <- colSums(sequence_PA)
samples_per_locus <- data.frame(locus = names(tmp), n_samples = tmp, check.names = F, row.names = NULL)
write.csv(samples_per_locus, paste0(path_to_out, "samples_per_locus.csv"))

# Calculate percentage of length recovered for each exon in each sample --------
percent_length = sweep(seq_lengths, 2, reference_length, "/")
percent_length <- apply(percent_length, c(1,2), FUN = function(x) ifelse(x>1,1,x)) # checks no value is above 1

# Prepare matrix for heatmap ----------------------------------------------
#set thresholds
limits <- seq(0.01, 0.99, 0.01)

#make empty matrix
loci_stats <- matrix(nrow = length(limits), ncol = length(limits))

# loop through threshold combinations
pb <- txtProgressBar(min = 0, max = length(limits), style = 3, width = 100, char = "=")  
for (i in 1:length(limits)) {
  for (j in 1:length(limits)) {
    
    #reset limits
    percent_len_limit <- percent_length
    
    # if percentage exon length recovered is >= limit make value 1
    # if not make value 0
    percent_len_limit[which(percent_len_limit >= limits[i])] <-  1
    percent_len_limit[which(percent_len_limit < limits[i])] <-  0
    
    # For each locus calculate number of samples with >= limit
    col <- colSums(percent_len_limit != 0)
    
    # Calculate % samples with >= limits of exon for each exon
    col <- col/nrow(percent_length)
    
    loci_stats[i, j] <-  table(col >= limits[j])["TRUE"]
    
    setTxtProgressBar(pb, i)
    # print(paste(i, "% of exon length in ", j, " % of samples", sep=""))
  }
}
close(pb)

# Gradient plot (heatmap) -------------------------------------------------
loci_stats <- data.frame(loci_stats)

for (i in 1:length(limits)) {
  colnames(loci_stats)[i] <-  limits[i] * 100
  loci_stats$percent_of_locus_length[i] <-  limits[i] * 100
}

#melt table
loci_stats_melt <- melt(loci_stats, id.vars = "percent_of_locus_length")
# head(loci_stats_melt)
loci_stats_melt$variable <- as.numeric(loci_stats_melt$variable)
colnames(loci_stats_melt) <- c("percent_of_locus_length", "percent_of_samples", "n_loci")

# Gradient color
pal <- wes_palette("Zissou1", 100, type = "continuous")
p <- ggplot(loci_stats_melt, aes(percent_of_locus_length, percent_of_samples, z = n_loci))
p + geom_raster(aes(fill = n_loci)) +
  scale_fill_gradientn(colours = pal, name = "No. exons", n.breaks = 10) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(y= "% samples with exon", x = "% exon length recovered", caption = paste0("Max. ", max(loci_stats_melt$n_loci, na.rm = T), " exons were recovered (over the ", length(reference_length), " in the target sequences set).")) +
  geom_contour(colour = "white", linetype = "dashed", alpha = 0.5, bins = 6) +
  annotate(geom = "point", x = 100*limit_perc_length_wanted, y = 100*limit_perc_nb_wanted, color = "grey40")+
  annotate(geom = "text", x = 100*limit_perc_length_wanted, y = 100*limit_perc_nb_wanted, label = loci_stats_melt$n_loci[which(loci_stats_melt$percent_of_locus_length == (limit_perc_length_wanted*100) & loci_stats_melt$percent_of_samples == (limit_perc_nb_wanted*100))], hjust =-0.5, color = "grey40")+
  theme(
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 11, margin = margin(
      t = 5,
      r = 0,
      b = 0,
      l = 0
    )),
    legend.key = element_blank(),
    legend.title = element_text(),
    legend.text = element_text(size = 10),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank()
  )

ggsave(paste0(path_to_out, "exon_recovery_gradient_samples.png"), width = 8, height = 8)


# Save number of locus recovered in a table -------------------------------

limit_perc_length_wanted = c(0.5, 0.5, 0.5, 0.5)
limit_perc_nb_wanted = c(0.1, 0.20, 0.5, 0.75)
#set thresholds
limits <- c(0,0.01, 0.1, 0.20, 0.25, 0.5, 0.75, 0.80, 0.9,0.99, 1)

#make empty matrix
loci_stats <- matrix(nrow = length(limits), ncol = length(limits))

#loop through threshold combinations
# = same loop as above
for (i in 1:length(limits)) {
  for (j in 1:length(limits)) {
    percent_len_limit <- percent_length
    percent_len_limit[which(percent_len_limit >= limits[i])] <-  1
    percent_len_limit[which(percent_len_limit < limits[i])] <-  0
    col <- colSums(percent_len_limit != 0)
    col <- col/nrow(percent_length)
    loci_stats[i, j] <-  table(col >= limits[j])["TRUE"]
  }
}

#make table presentable
colnames(loci_stats) <- paste0(limits*100, "% samples")
rownames(loci_stats) <- paste0(limits*100, "% exons")
loci_stats[is.na(loci_stats)] <- 0
write.csv(loci_stats, paste0(path_to_out, "exon_filtering_stats.csv"))


# Filter the list of loci according to desired percentage -----------------
# NEEDS TO DEFINE THE VALUES OF limit_perc_length_wanted AND limit_perc_nb_wanted BEFOREHAND
# NEEDS TO HAVE limit_perc_length_wanted AND limit_perc_nb_wanted TO BE THE SAME LENGTH
# otherwise default will be set 0.75 for both
if (! "limit_perc_length_wanted" %in% ls() || length(limit_perc_length_wanted) != length(limit_perc_nb_wanted)) limit_perc_length_wanted <- 0.50
if (! "limit_perc_nb_wanted" %in% ls() || length(limit_perc_length_wanted) != length(limit_perc_nb_wanted)) limit_perc_nb_wanted <- 0.20

for (i in 1:length(limit_perc_length_wanted)){
  # if percentage exon length recovered is >= limit_perc_length_wanted make value 1
  # if not make value 0
  percent_len_limit <- percent_length
  percent_len_limit[which(percent_len_limit >= limit_perc_length_wanted[i])] = 1
  percent_len_limit[which(percent_len_limit < limit_perc_length_wanted[i])] = 0
  
  # For each exon, calculate number of samples in which the exon is recovered with a length >= limit_perc_length_wanted
  col <- colSums(percent_len_limit)
  
  # From this number, calculate the % of samples (in which the exon is recovered with a length >= limit_perc_length_wanted)
  col <- col/(nrow(seq_lengths)) # percentage 
  # table(col >= limit_perc_nb_wanted[i])
  
  # Retrieve the list of exons that are recovered with a length greater than limit_perc_length_wanted in more than limit_perc_nb_wanted samples
  filtered_exons <- names(which(col >= limit_perc_nb_wanted[i]))
  
  # Filter the original dataset (exon length) to keep only the filtered exons
  length_wanted <- seq_lengths[, filtered_exons]
  
  # Density plot of filtered exon lengths -----------------------------------
  #Density plot of exon lengths with rug of actual exon lengths
  #make data frame
  marker_len <- data.frame(filtered_exons, as.numeric(reference_length[filtered_exons]))
  colnames(marker_len) <- c("locus", "length")
  
  # Basic density plot with mean line and marginal rug
  ggdensity(
    marker_len,
    x = "length",
    fill = "#0073C2FF",
    color = "#0073C2FF",
    alpha = 0.15,
    add = "mean",
    rug = TRUE,
    xlab = "Recovered exon length (bp)",
    ylab = "Density",
    title = paste0(100*limit_perc_length_wanted[i],"_",100*limit_perc_nb_wanted[i], " filter"),
    xlim = c(0, max(reference_length))
  )
  
  ggsave(paste0(path_to_out, (100*limit_perc_length_wanted[i]), "_", (100*limit_perc_nb_wanted[i]), "_exon_length_density.png"))
  
  
  # Create copy commands for keep filtered loci - for use downstream --------
  cat( 
    x = paste0("mv *", filtered_exons, "*.* ", "filtered","/"),
    file = paste0(path_to_out, "move_", (100*limit_perc_length_wanted[i]), "_", (100*limit_perc_nb_wanted[i]), ".sh"),
    sep = "\n"
  )
  
  # Create list of filtered loci - for use downstream --------
  cat(
    x = sort(paste0(filtered_exons)),
    file = paste0(path_to_out, "list_", (100*limit_perc_length_wanted[i]), "_", (100*limit_perc_nb_wanted[i]), ".txt"),
    sep = "\n"
  )
}

