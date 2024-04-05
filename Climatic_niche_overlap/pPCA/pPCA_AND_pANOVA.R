#########################################################################################################
### Phylogenetic PCA and Phylogenetic ANOVA. Scripts modified from Liam Revell Phytools blog http://blog.phytools.org/search?q=phyl.pca
#########################################################################################################

require("rJava")
library("rgdal")
library("spatial")
library("tidyverse")
library("dismo")
library("raster")
library("ENMeval")
library("stringr")
library("sp")
library("GGally")
library("ggplot2")
library("sf")
library("BIEN")
library("rgbif")
library("maptools")
library("dplyr")
library("rlist")
library("devtools")
library("ggfortify")
library("nicheROVER")
library("phytools")
library("ape")
library("ComplexHeatmap")
library("geiger")
library("phangorn")
library("phylogram")
library("circlize")
library('RColorBrewer')
library("factoextra")

setwd("Climatic_niche_overlap/pPCA/")
#source CANDI functions
source("~/Applications/candi/R/get_occ_records.R")
source("~/Applications/candi/R/get_both_occ_records.R")
source("~/Applications/candi/R/get_bien_occ_records.R")
source("~/Applications/candi/R/get_gbif_occ_records.R")
source("~/Applications/candi/R/remove_dup_locs.R")
source("~/Applications/candi/R/remove_ocean_points.R")
source("~/Applications/candi/R/remove_perf_0_90_180.R")
source("~/Applications/candi/R/remove_lessthan.R")
source("~/Applications/candi/R/remove_null_items.R")
source("~/Applications/candi/R/get_world_clim.R")
source("~/Applications/candi/R/make_corr_matrix.R")
source("~/Applications/candi/R/remove_corr_variables.R")
source("~/Applications/candi/R/remove_points_outside_nat_range.R")
source("~/Applications/candi/R/model_niches.R")
source("~/Applications/candi/R/get_world_clim.R")
source("~/Applications/candi/R/trim_to_shapefile.R")

## CLIMATE DATA PREP PATHWAY ##
###############################

load("~/Palicourea/Climatic_niche_overlap/occurrences_filtered.RData")

coords<-c("longitude","latitude")
coordinates = lapply(occurrences, "[", , coords)
list2env(coordinates,envir=.GlobalEnv)

#download 19 bioclim variables from WorldClim and process into a RasterStack
climate_stack <- get_world_clim()

#makes pdf of correlation matrix saved to working directory to pick out correlated variables
#matrix <- make_corr_matrix(occurrences = occurrences, environment_data = climate_stack)

#select variables to remove
bad_vars <- c("bio10", "bio11", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18",
              "bio19", "bio2", "bio3", "bio4", "bio7", "bio8", "bio9")

#remove variables selected by user
uncorr_stack <- remove_corr_variables(raster_stack = climate_stack, variables_to_be_removed = bad_vars) 
final_stack <- uncorr_stack

##EXTRACTING BIOCLIM VALUES FROM OCCURRENCE POINTS##
####################################################
points <- lapply(coordinates, function(x) SpatialPoints(x, proj4string = final_stack@crs))
values <- lapply(points, function(x) extract(final_stack,x))

#bind lists into df
values_df<-list.rbind(values)
log_values_df<-log10(values_df)
spp_df<-list.rbind(occurrences)

df <- cbind.data.frame(spp_df,log_values_df)

niche_data<-df[4:8]
niche_data$species<-gsub(" ", "_", niche_data$species)
niche_data_complete<-na.omit(niche_data)

##READ TREE #########################################
#####################################################

tree<- ape::read.tree("~/Palicourea/Climatic_niche_overlap/pPCA/dendrogram.tre")

tips2delete<-c("Car_guianensis_Gonzalez2158","Carapichea_ipecacuanha_Croat15117",
               "Eum_boliviana_Campbell22035","Not_epiphytica_Neill15737",
               "Not_uliginosa_Stevens37138","Psy_carthagenensis_Araujo2124",
               "Psy_grandis_Taylor11745","Psy_guianensis_Merello1711",
               "Psy_horizontalis_Stevens32733","Psy_jinotegensis_Stevens33549",
               "Psy_limonensis_Stevens31580","Psy_marginata_Stevens32781",
               "Psy_nervosa_Stevens32362","Psy_panamensis_Stevens32285","Psy_subsessilis_Stevens31494","Rud_cornifolia_deGracias818")

tree <- ape::drop.tip(tree, tips2delete)

group_mean<- aggregate(x= niche_data_complete[2:5],
                       # Specify group indicator
                       by = list(niche_data_complete$species),      
                       # Specify function (i.e. mean)
                       FUN = mean)
areas<-c("L","E","A","A","E","A","A","A","A","L","F","M","E","A",
         "C","E","E","E","C","E","A","E","E","F","A","E","C","A",
         "A","E","L","E","E","E","C","A","E","E","A","A","A","A",
         "E","E","E","E","E","A","A","E","E","A","A","E","M","A",
         "E","E","A","F","A","A","A","M","E","A","A","A","L","A",
         "L","M","E","E","E","L","F") ##clumsy but works. You make it better
group_mean$areas<-areas


#turn species' column into names
group_mean_names<-row.names(group_mean) <- group_mean$Group.1
group_mean_names[1] <- NULL

###Phylogenetic PCA##############################################################
#################################################################################

pPCA<-phyl.pca(tree,group_mean_names[1:4])
obj<-as.princomp(pPCA)

fviz_screeplot(obj,addlabels=TRUE)
fviz_pca_ind(obj,label="none",habillage=areas,
             addEllipses=F)


###Phylogenetic ANOVA############################################################
#################################################################################
#drop values for species outside A and E

group_meanAE<-group_mean[(group_mean$area =="A")| (group_mean$area=="E"),]
species<-group_meanAE$Group.1

#extract values for each variable and add labels (nneds loop)
bio1<-group_meanAE$bio1
bio1<- setNames(bio1, species)
areas<-group_meanAE$areas
areas<-setNames(areas, species)
bio12<-group_meanAE$bio12
bio12<- setNames(bio12, species)
bio5<-group_meanAE$bio5
bio5<- setNames(bio5, species)
bio6<-group_meanAE$bio6
bio6<- setNames(bio6, species)

##remove tips that are not in A or E
tips2delete<-c("Palicourea_acanthacea","Palicourea_brachiata","Palicourea_brachypoda",
               "Palicourea_brevicollis","Palicourea_correae","Palicourea_correae","Palicourea_cyanococca",
               "Palicourea_divaricata","Palicourea_elata","Palicourea_glomerulata","Palicourea_hazenii",
               "Palicourea_racemosa","Palicourea_sessilis","Palicourea_subfusca","Palicourea_timbiquensis",
               "Palicourea_tomentosa","Palicourea_topoensis","Psychotria_rosea","Psychotria_suterella")
tree <- ape::drop.tip(tree, tips2delete)

#phylANOVAs
phylANOVA(tree,areas,bio1)
phylANOVA(tree,areas,bio12)
phylANOVA(tree,areas,bio5)
phylANOVA(tree,areas,bio6)

