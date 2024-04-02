library("phytools")
library("phyloch")
library("strap")
library("coda")

setwd("~/Desktop/Palicourea/RevBayes/simple/")
tips2delete<-c("Car_guianensis_Gonzalez2158","Carapichea_ipecacuanha_Croat15117",
               "Eum_boliviana_Campbell22035","Not_epiphytica_Neill15737",
               "Not_uliginosa_Stevens37138","Psy_carthagenensis_Araujo2124",
               "Psy_grandis_Taylor11745","Psy_guianensis_Merello1711",
               "Psy_horizontalis_Stevens32733","Psy_jinotegensis_Stevens33549",
               "Psy_limonensis_Stevens31580","Psy_marginata_Stevens32781",
               "Psy_nervosa_Stevens32362","Psy_panamensis_Stevens32285","Psy_subsessilis_Stevens31494","Rud_cornifolia_deGracias818")

#mcc<-"~/Desktop/Palicourea/Manuscript/calibration/standard_Palicourea_ORC_versioncorrected_MCC_newick.tre"

tree<-read.beast("~/Desktop/Palicourea/Manuscript/calibration/standard_Palicourea_ORC_versioncorrected_MCC.tre")

t <- ape::drop.tip(tree, tips2delete)

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


pdf("tree.pdf", width = 20, height = 20)

geoscalePhylo(tree = t,
              x.lim = c(-2,21),
              units = c("Epoch"),
              tick.scale = "myr",
              boxes = FALSE,
              width = 1,
              cex.tip = 2,
              cex.age = 3,
              cex.ts = 2,
              erotate = 0,
              label.offset = 0.1)

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

dev.off()

plotTree(tree,xlim=c(110,-5),direction="leftwards",
         mar=c(4.1,1.1,1.1,1.1),ftype="i")
abline(v=seq(0,120,by=10),lty="dashed",
       col=make.transparent("grey",0.5))
axis(1,at=seq(0,120,by=20))
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
for(i in 1:tree$Nnode+Ntip(tree))
  lines(x=c(CI[i-Ntip(tree),1],CI[i-Ntip(tree),2]),
        y=rep(obj$yy[i],2),lwd=11,lend=0,
        col=make.transparent("blue",0.4))
points(obj$xx[1:tree$Nnode+Ntip(tree)],
       obj$yy[1:tree$Nnode+Ntip(tree)],pch=19,col="blue",
       cex=1.8)
#write.tree(tree,file="rev_dendrogram.tre")
