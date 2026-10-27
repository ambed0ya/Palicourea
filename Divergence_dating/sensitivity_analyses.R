library(RevGadgets)
library(ape)
library(phangorn)
library(ggplot2)
library(dplyr)

# MCC tree used to define node 227
mcc_file <- "~/Desktop/Palicourea_temp/prioronly_mcc.tre"

posterior_files <- c(
  "~/Desktop/Palicourea_temp/UCLN_prioronly_combined.trees",
  "~/Desktop/Palicourea_temp/UCLN_yule_combined.trees",
  "~/Desktop/Palicourea_temp/UCLN_combined.trees"
)

labels <- c("Only prior","Yule", "Posterior")


# get clade tips for Palicoureeae
mcc_obj <- readTrees(paths = mcc_file)
mcc_phy <- mcc_obj[[1]][[1]]@phylo

target_node <- getMRCA(mcc_phy,c("Pal_stipulosa", "Pal_lehmannii"))

desc <- phangorn::Descendants(mcc_phy, target_node, type = "tips")[[1]]
clade_tips <- mcc_phy$tip.label[desc]

length(clade_tips)
clade_tips

all_dfs <- list()

for (i in seq_along(posterior_files)) {
  
  cat("\nProcessing", i, "of", length(posterior_files),
      "-", labels[i], "\n")
  
  post_obj <- readTrees(paths = posterior_files[i])
  post_trees <- post_obj[[1]]
  
  cat("Trees:", length(post_trees), "\n")
  
  # These are the descendants of MCC node 227 (Palicourea MRCA)
  # that occur in the posterior trees
  phy1 <- post_trees[[1]]@phylo
  tips_present <- clade_tips[clade_tips %in% phy1$tip.label]
  
  cat("Target tips:", length(tips_present),
      "of", length(clade_tips), "\n")
  
  heights <- vapply(post_trees, function(tree_obj) {
    
    phy <- tree_obj@phylo
    
    # Find the MRCA of the node-227 descendants
    # in THIS posterior tree
    node <- ape::getMRCA(phy,tips_present)
    
    if (is.null(node) || is.na(node))
      return(NA_real_)
    
    # Distance from root to each node
    depth <- ape::node.depth.edgelength(phy)
    
    # Total height of tree
    tree_height <- max(
      depth[seq_len(ape::Ntip(phy))]
    )
    
    # Age of the MRCA corresponding to MCC node 227
    tree_height - depth[node]
    
  }, numeric(1))
  
  all_dfs[[i]] <- data.frame(
    height = heights,
    analysis = labels[i]
  ) %>%
    filter(!is.na(height))
  
  cat("Done:", labels[i], "\n")
}

plot_df <- bind_rows(all_dfs)


##Plotting results
ggplot(plot_df, aes(x = height, fill = analysis, color = analysis, alpha = analysis)) +
  geom_density(linewidth = 0.8) +
  scale_alpha_manual(values = c(
    "Only prior" = 0.35,
    "Yule" = 0.20,
    "Posterior" = 0.35
  )) +
  theme_classic() +
  labs(
    x = "Node age",
    y = "Density",
    fill = "Analysis",
    color = "Analysis"
  )


