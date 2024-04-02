Modified from David Černý (https://davidcerny.github.io/post/plotting_beast/)

library("phytools")
library("phyloch")
library("strap")
library("coda")

annot_tree <- phyloch::read.beast("~/Desktop/Palicourea/Manuscript/calibration/standard_Palicourea_ORC_versioncorrected_MCC.tre")
#t<-read.beast("~/Desktop/Palicourea/Manuscript/calibration/standard_Palicourea_ORC_versioncorrected_MCC.tre")
annot_tree <- ape::drop.tip(annot_tree, tips2delete)

if (is.null(annot_tree$`CAheight_95%_HPD_MIN`)) {
  annot_tree$min_ages <- annot_tree$`height_95%_HPD_MIN`
  annot_tree$max_ages <- annot_tree$`height_95%_HPD_MAX`
} else {
  annot_tree$min_ages <- annot_tree$`CAheight_95%_HPD_MIN`
  annot_tree$max_ages <- annot_tree$`CAheight_95%_HPD_MAX`
}

annot_tree$root.time <- max(nodeHeights(annot_tree)) + 0.0

pdf("tree.pdf", width = 20, height = 20)
geoscalePhylo(ladderize(annot_tree, right = F), x.lim = c(-10, 45), cex.tip = 0.5, cex.age = 1.3, cex.ts = 0.4)

T1 <- get("last_plot.phylo", envir = .PlotPhyloEnv)

for(i in (Ntip(annot_tree) + 1):(annot_tree$Nnode + Ntip(annot_tree))) {
  lines(x = c(T1$root.time - annot_tree$min_ages[i - Ntip(annot_tree)],
              T1$root.time - annot_tree$max_ages[i - Ntip(annot_tree)]),
        y = rep(T1$yy[i], 2), lwd = 4, lend = 0,
        col = make.transparent("blue", 0.4))
}
dev.off()
setwd("~/Desktop/")
t$root.time <- t$height[1]

num_taxa <- length(t$tip.label)

display_all_node_bars <- TRUE

names_list <-vector()
for (name in t$tip){
  v <- strsplit(name, "_")[[1]]
  if(display_all_node_bars){
    names_list = c(names_list, name)
  }
  else if(v[length(v)]=="0"){
    names_list = c(names_list, name)
  }
}

nids <- vector()
pos <- 1
len_nl <- length(names_list)
for(n in names_list){
  for(nn in names_list[pos:len_nl]){
    if(n != nn){
      m <- getMRCA(t,c(n, nn))
      if(m %in% nids == FALSE){
        nids <- c(nids, m)
      }
    }
  }
  pos <- pos+1
}




geoscalePhylo(tree = t,
              x.lim = c(0,45),
              units = c("Epoch"),
              tick.scale = "myr",
              boxes = FALSE,
              width = 1,
              cex.tip = 2,
              cex.age = 3,
              cex.ts = 2,
              erotate = 0,
              label.offset = 0.01)

lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

for(nv in nids){
  bar_xx_a <- c(lastPP$xx[nv]+t$height[nv-num_taxa]-t$"height_95%_HPD_MIN"[nv-num_taxa],
                lastPP$xx[nv]-(t$"height_95%_HPD_MAX"[nv-num_taxa]-t$height[nv-num_taxa]))
  lines(bar_xx_a, c(lastPP$yy[nv], lastPP$yy[nv]), col = rgb(0, 0, 1, alpha = 0.3), lwd = 12)
}

t$node.label <- t$posterior
p <- character(length(t$node.label))
p[t$node.label >= 0.95] <- "black"
p[t$node.label < 0.95 & t$node.label >= 0.75] <- "gray"
p[t$node.label < 0.75] <- "white"
nodelabels(pch = 21, cex = 1.5, bg = p)



#plotTree(t,xlim=c(50,-5),direction="leftwards",
#         mar=c(4.1,1.1,1.1,1.1),ftype="i")
#abline(v=seq(0,50,by=5),lty="dashed",
#       col=make.transparent("grey",0.5))
#axis(1,at=seq(0,120,by=20))
#obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
#for(i in 1:t$Nnode+Ntip(t))
#  lines(x=c(CI[i-Ntip(t),1],CI[i-Ntip(t),2]),
#        y=rep(obj$yy[i],2),lwd=11,lend=0,
#        col=make.transparent("blue",0.4))
#points(obj$xx[1:t$Nnode+Ntip(tree)],
#       obj$yy[1:t$Nnode+Ntip(tree)],pch=19,col="blue",
#       cex=1.8)
#write.tree(tree,file="rev_dendrogram.tre")
