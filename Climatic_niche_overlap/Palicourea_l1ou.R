# Palicourea L1OU analysis
setwd("~/Desktop/Submitted/Palicourea/Manuscript/niche/")
library(devtools)
library(l1ou)
library(ape)
library(treeio)
library(phangorn)
library(phytools)

astral01 <- read.tree("rev_dendrogram.tre")
tips2delete<-c("Pal_acuminata1","Pal_guianensis2","Pal_quinquepyrena")
tree <- ladderize(ape::drop.tip(astral01, tips2delete), right = F)
plot(tree)

# Update tip labels so that they match the names in the niche data
new_names <- read.csv("pal_name_updates.csv", header = T)
match_indices <- match(tree$tip.label, new_names$old_label)
# Replace names in tree$tip.label vector
tree$tip.label <- ifelse(is.na(match_indices), tree$tip.label, new_names$new_label[match_indices])
plot(tree)

niche_dat <- read.table("Palicourea_median_climate_data_with_elev.csv", header = T, sep = "\t", row.names = "species")[,2:6] #all vars separate
head(niche_dat, 10)

pal <- adjust_data(tree, as.matrix(niche_dat))
colnames(pal$Y) <- c("bio1", "bio12", "bio5", "bio6", "elev")
eModel <- estimate_shift_configuration(pal$tree, pal$Y, criterion = "pBIC")

# Set branch/shift colors
# Plot phylo with edge labels for reference
plot(eModel$tree)
edgelabels(seq_along(tree$edge[, 1]), adj = c(0.5, -0.5), frame = "n", cex = 0.8, col = "blue")
# Make new df correlating shift configurations with range/color with corresponding node at cladogenesis 
shift_cols <- data.frame(node = eModel$shift.configuration, 
                         colors = c("deepskyblue", "#568259", "deepskyblue", "deepskyblue", "#97CC04", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue",
                                    "#FB9B2D", "yellow", "deepskyblue", "deepskyblue", "yellow", "deepskyblue", "deepskyblue", "deepskyblue", "yellow", "#568259"))

pdf(file = "l1ou.pdf", width = 16, height = 12)
plot(eModel, palette = c(shift_cols$colors, "gray90"), 
     asterisk = T, edge.shift.ann = F, 
     cex = 1.5, edge.width = 3, label.offset = 0.05, 
     bar.axis = T, show.tip.label = F)
dev.off()

#result <- l1ou_bootstrap_support(eModel, nItrs=100, multicore=TRUE, nCores=4)
#result$detection.rate
#plot(result)
#result$detection.rate


