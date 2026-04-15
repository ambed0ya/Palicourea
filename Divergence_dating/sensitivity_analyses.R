#NODE POSTERIORS FOR ALL TREES
library(RevGadgets)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(viridis)

#FORCED ROOT
# Define base path
setwd("~/Downloads/")
base_path <- "~/Downloads/"

# Define tree paths
tree_paths <- c(
  "HKY_onlyprior_mcc.tree",
  "HKY_sd6_mcc.tre",
  "HKY_yule_mcc.tre",
  "HKY_fixed_dedupIDs_strict_MCC.tre"
)

# Descriptive plot titles
plot_titles <- c(
  "Only prior",
  "SD6",
  "Yule",
  "posterior"
)

# Empty list to store plots
all_plots <- list()

# Loop through files and make plots
for (i in seq_along(tree_paths)) {
  full_path <- file.path(base_path, tree_paths[i])
  tree <- readTrees(paths = full_path)
  tree_data <- tree[[1]][[1]]
  annotations <- tree_data@data
  
  expanded <- annotations %>%
    filter(!is.na(age_0.95_HPD)) %>%
    rowwise() %>%
    mutate(posterior_values = list(seq(age_0.95_HPD[1], age_0.95_HPD[2], length.out = 100))) %>%
    unnest(posterior_values)
  
  p <- ggplot(expanded, aes(x = posterior_values, fill = as.factor(node))) +
    geom_density(alpha = 0.5, adjust = 0.5) +
    labs(title = plot_titles[i], x = NULL, y = NULL) +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(size = 10, hjust = 0.5), panel.grid.minor = element_blank()) +
    scale_fill_viridis_d() +
    scale_x_reverse(limits = c(11, -1), breaks = 0:10) + 
    coord_cartesian(ylim = c(0, 2.25))
  
  all_plots[[i]] <- p
}

# Reordered plot list
ordered_plots <- c(
  all_plots[1], all_plots[6], all_plots[2], all_plots[7], all_plots[3], # Row 1 (no migration)
  all_plots[8], all_plots[4], all_plots[9], all_plots[5], all_plots[10]  # Row 2 (with gene flow)
)

# Arrange all plots in a grid (2 columns x 5 rows)
grid.arrange(grobs = ordered_plots, ncol = 2)
