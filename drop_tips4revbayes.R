library(geiger)
library(ape)
library(phytools)
library(phangorn)
library(TreeTools)
library(phylogram)

setwd("~/Desktop/Palicourea/RevBayes/simple/")
tips2delete<-c("Car_guianensis_Gonzalez2158","Carapichea_ipecacuanha_Croat15117",
               "Eum_boliviana_Campbell22035","Not_epiphytica_Neill15737",
               "Not_uliginosa_Stevens37138","Psy_carthagenensis_Araujo2124",
               "Psy_grandis_Taylor11745","Psy_guianensis_Merello1711",
               "Psy_horizontalis_Stevens32733","Psy_jinotegensis_Stevens33549",
               "Psy_limonensis_Stevens31580","Psy_marginata_Stevens32781",
               "Psy_nervosa_Stevens32362","Psy_panamensis_Stevens32285","Psy_subsessilis_Stevens31494","Rud_cornifolia_deGracias818")

trees<-"~/Desktop/Palicourea/Manuscript/calibration/standard_Palicourea_ORC_versioncorrected_MCC_newick.tre"
trees
trees01<-ape::read.tree(trees)
trees01
tree <- ape::drop.tip(trees01, tips2delete)
plot(tree)
write.tree(tree,file="rev_dendrogram.tre")
