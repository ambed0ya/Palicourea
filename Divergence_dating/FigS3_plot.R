Modified from David Černý (https://davidcerny.github.io/post/plotting_beast/)

library("phytools")
library("phyloch")
library("strap")
library("coda")
setwd("~/Palicourea/Divergence_dating/")
annot_tree <- phyloch::read.beast("standard_Palicourea_ORC_MCC.tre")
#t<-read.beast("~/Desktop/Palicourea/Manuscript/calibration/standard_Palicourea_ORC_versioncorrected_MCC.tre")
#annot_tree <- ape::drop.tip(annot_tree, tips2delete)

if (is.null(annot_tree$`CAheight_95%_HPD_MIN`)) {
  annot_tree$min_ages <- annot_tree$`height_95%_HPD_MIN`
  annot_tree$max_ages <- annot_tree$`height_95%_HPD_MAX`
} else {
  annot_tree$min_ages <- annot_tree$`CAheight_95%_HPD_MIN`
  annot_tree$max_ages <- annot_tree$`CAheight_95%_HPD_MAX`
}

annot_tree$root.time <- max(nodeHeights(annot_tree)) + 0.0

pdf("Fig.S3.pdf", width = 100, height = 100)
geoscalePhylo(ladderize(annot_tree, right = F), x.lim = c(-10, 45), cex.tip = 6, cex.age = 10, cex.ts = 9)

T1 <- get("last_plot.phylo", envir = .PlotPhyloEnv)

for(i in (Ntip(annot_tree) + 1):(annot_tree$Nnode + Ntip(annot_tree))) {
  lines(x = c(T1$root.time - annot_tree$min_ages[i - Ntip(annot_tree)],
              T1$root.time - annot_tree$max_ages[i - Ntip(annot_tree)]),
        y = rep(T1$yy[i], 2), lwd = 30, lend = 0,
        col = make.transparent("blue", 0.4))
}
dev.off()
