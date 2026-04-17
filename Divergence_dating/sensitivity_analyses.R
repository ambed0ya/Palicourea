library(RevGadgets)
library(ape)
library(phangorn)
library(ggplot2)
library(dplyr)

# MCC tree used to define node 227
mcc_file <- "~/Downloads/HKY_onlyprior_mcc.tree"

posterior_files <- c(
  "~/Downloads/HKY_dedupIDs_onlyprior.trees",
  "~/Downloads/HKY_dedupIDs_strict_yule.trees",
  "~/Downloads/HKY_fixed_dedupIDs_strict_combined.trees"
)

labels <- c("Only prior","Yule", "Posterior")


# get clade tips for Palicoureeae
mcc_obj <- readTrees(paths = mcc_file)
mcc_phy <- mcc_obj[[1]][[1]]@phylo

getMRCA(mcc_phy, c("Psy_stipulosa_Zappi984", "Pal_lehmannii_Silverstone8396"))
target_node <- 227
desc <- phangorn::Descendants(mcc_phy, target_node, type = "tips")[[1]]
clade_tips <- mcc_phy$tip.label[desc]

get_clade_height <- function(tree_obj, tips) {
  phy <- tree_obj@phylo
  tips_present <- intersect(tips, phy$tip.label)
  
  if (length(tips_present) < 2) return(NA_real_)
  
  node <- ape::getMRCA(phy, tips_present)
  if (is.null(node) || is.na(node)) return(NA_real_)
  
  bt <- ape::branching.times(phy)
  if (!(as.character(node) %in% names(bt))) return(NA_real_)
  
  as.numeric(bt[as.character(node)])
}

all_dfs <- list()

for (i in seq_along(posterior_files)) {
  post_obj <- readTrees(paths = posterior_files[i])
  post_trees <- post_obj[[1]]
  
  heights <- sapply(post_trees, get_clade_height, tips = clade_tips)
  
  all_dfs[[i]] <- data.frame(
    height = heights,
    analysis = labels[i]
  ) %>%
    filter(!is.na(height))
}

plot_df <- bind_rows(all_dfs)


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


