#########################################################################################################
### 1.Extracting and cleaning occurrence and environmental data for large datasets with CANDI pipeline###
# Modified from Tribble et al. 2023.
#########################################################################################################
setwd("~/Desktop/Submitted/Palicourea/Manuscript/niche")

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
library('factoextra')

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

## OCCURRENCE DATA PREP PATHWAY ##
##################################

#load data# 
#world map
big_map <- rgdal::readOGR("gadm28_adm0/gadm28_adm0.shp")
#species range data
read_csv("Palicourea_native_ranges.csv") %>%
  dplyr::select(species, range) -> all_native_ranges
#botanical map not used here 
#kew_map_level_2 <- rgdal::readOGR("wgsrpd-master/level2/level2.shp")

#species list 
my_species1<-c("Palicourea acanthacea","Palicourea acuminata","Palicourea allenii","Palicourea amethystina","Palicourea andina")
my_species2<-c("Palicourea angustifolia","Palicourea apicata","Palicourea bangii","Palicourea brachiata","Palicourea brevicollis","Palicourea lasiantha")
my_species3<-c("Palicourea conephoroides","Palicourea correae","Palicourea corymbifera","Palicourea crocea","Palicourea croceoides", "Palicourea cyanococca","Palicourea demissa","Palicourea didymocarpos")
my_species4<-c("Palicourea divaricata","Palicourea domingensis","Palicourea egensis","Palicourea flavescens","Palicourea flavifolia","Palicourea glomerulata")
my_species5<-c("Palicourea grandiflora", "Palicourea guianensis","Palicourea hazenii","Palicourea lehmannii","Palicourea lineata","Palicourea loxensis","Palicourea luteonivea")
my_species6<-c("Palicourea macrobotrys","Palicourea marcgravii","Palicourea padifolia","Palicourea petiolaris","Palicourea polycephala","Palicourea pyramidalis","Palicourea quadrifolia")
my_species7<-c("Palicourea quinquepyrena","Palicourea racemosa","Palicourea rhodothamna","Palicourea rigida","Palicourea seemannii","Palicourea standleyana","Palicourea stenosepala")
my_species8<-c("Palicourea stipularis","Palicourea subfusca","Palicourea sulphurea","Palicourea tetragona","Palicourea thyrsiflora","Palicourea timbiquensis","Palicourea tinctoria")
my_species9<-c("Palicourea topoensis","Palicourea triphylla","Palicourea sessilis","Palicourea woronovii","Palicourea brachypoda","Palicourea berteroana","Palicourea winkleri")
my_species10<-c("Palicourea callithrix","Palicourea deflexa","Palicourea elata","Palicourea gracilenta","Palicourea justiciifolia","Palicourea ostreophora","Palicourea obliquinervia", "Psychotria rosea")
my_species11<-c("Palicourea dichotoma","Palicourea tomentosa","Palicourea prunifolia","Palicourea pubescens","Palicourea reticulata","Palicourea suerrensis")
my_species12<-c("Psychotria suterella","Palicourea glabra","Palicourea_jelskii","Palicourea_nitidella")



#download occurrence data from "bien" AND "gbif".
occ_data1 <- get_occ_records(species = my_species1, database = "both")
occ_data2 <- get_occ_records(species = my_species2, database = "both")
occ_data3 <- get_occ_records(species = my_species3, database = "both")
occ_data4 <- get_occ_records(species = my_species4, database = "both")
occ_data5 <- get_occ_records(species = my_species5, database = "both")
occ_data6 <- get_occ_records(species = my_species6, database = "both")
occ_data7 <- get_occ_records(species = my_species7, database = "both")
occ_data8 <- get_occ_records(species = my_species8, database = "both")
occ_data9 <- get_occ_records(species = my_species9, database = "both")
occ_data10 <- get_occ_records(species = my_species10, database = "both")
occ_data11 <- get_occ_records(species = my_species11, database = "both")
occ_data12 <- get_occ_records(species = my_species12, database = "both")

occurrences<-c(occ_data1,occ_data2,occ_data3,occ_data4,occ_data5,occ_data6,occ_data7,
               occ_data8,occ_data9,occ_data10,occ_data11,occ_data12)
save(occurrences, file = "~/Desktop/Palicourea/Manuscript/niche/occurrences_raw.RData")

#data cleaning
occurrences <- remove_dup_locs(occurrences)
occurrences <- remove_ocean_points(occurrences, world_map = big_map)
occurrences <- remove_perf_0_90_180(occurrences)
#occ_data <- remove_points_outside_nat_range(df = occurrences, 
#                                           botan_map = kew_map_level_2, 
#                                           nat_range = all_native_ranges)
occurrences <- remove_lessthan(occurrences, n = 10)
occurrences <- remove_null_items(occurrences)

save(occurrences, file = "~/Desktop/Palicourea/Manuscript/niche/occurrences_filtered.RData")

## CLIMATE DATA PREP PATHWAY ##
###############################

load("~/Desktop/Submitted/Palicourea/Manuscript/niche/occurrences_filtered.RData")

coords<-c("longitude","latitude")
coordinates = lapply(occurrences, "[", , coords)
list2env(coordinates,envir=.GlobalEnv)

#download 19 bioclim variables from WorldClim and process into a RasterStack
climate_stack <- get_world_clim()

#If the above does not work then download manually

#makes pdf of correlation matrix saved to working directory to pick out correlated variables
matrix <- make_corr_matrix(occurrences = occurrences, environment_data = climate_stack)

#select variables to remove
bad_vars <- c("bio10", "bio11", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18",
                             "bio19", "bio2", "bio3", "bio4", "bio7", "bio8", "bio9")

#remove variables selected by user
uncorr_stack <- remove_corr_variables(raster_stack = climate_stack, variables_to_be_removed = bad_vars) 
final_stack <- uncorr_stack
#final_stack <- climate_stack


##EXTRACTING BIOCLIM VALUES FROM OCCURRENCE POINTS##
####################################################
points <- lapply(coordinates, function(x) SpatialPoints(x, proj4string = final_stack@crs))
values <- lapply(points, function(x) extract(final_stack,x))

#bind lists into df
values_df<-list.rbind(values)
log_values_df<-log10(values_df)
spp_df<-list.rbind(occurrences)

df <- cbind.data.frame(spp_df,log_values_df)

#add elevation using raster manually downloaded from worldclim.org
elev <- raster("elevation.tif")

long_lat <- df[c("longitude", "latitude")]
coordinates(long_lat) <- ~longitude + latitude
crs(long_lat) <- crs(elev)
df$elev <- extract(elev, long_lat)

#alt method; extracts elevation automatically, but takes a really long time to download data
#pal_sf <- sf::st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326)
#df$elev <- get_elev_point(locations = pal_sf, prj = 4326)

niche_data<-df[4:9]
niche_data$species<-gsub(" ", "_", niche_data$species)
niche_data_complete<-na.omit(niche_data)
write.csv(niche_data_complete, "Palicourea_climate_data_with_elev.csv", row.names = F)

##Calculationn of median value for cimatic variables for each species
median_niche_values<-aggregate(niche_data_complete[2:6], niche_data_complete[1], median)
write.csv(median_niche_values, "Palicourea_median_climate_data_with_elev.csv", row.names=F)
#median_niche_values$species<-as.factor(median_niche_values$species)
#load data
median_niche_values<-read.csv("Palicourea_median_climate_data_with_elev.csv")

#########################################################################################
##################DISPARITY ACROSS GROUPS AND THROUGH TIME###############################
#########################################################################################
#convert Species column into row names
median_niche_values_rownames <- data.frame(median_niche_values[,-1], row.names=median_niche_values[,1])

#Convert dataframe to matrix
matrix_median_niche_values<-data.matrix(median_niche_values_rownames)
#Calculate and plot DTT across variables
par( mfrow= c(5,1) )
dttbio1<-dtt(phy=tree, data=matrix_median_niche_values[,"bio1"], index="avg.sq", nsim=1000, plot=TRUE)
dttelev<-dtt(phy=tree, data=matrix_median_niche_values[,"elev"], index="avg.sq", nsim=1000, plot=TRUE)
dttbio6<-dtt(phy=tree, data=matrix_median_niche_values[,"bio6"], index="avg.sq", nsim=1000, plot=TRUE)
dttbio12<-dtt(phy=tree, data=matrix_median_niche_values[,"bio12"], index="avg.sq", nsim=1000, plot=TRUE)
dttbio5<-dtt(phy=tree, data=matrix_median_niche_values[,"bio5"], index="avg.sq", nsim=1000, plot=TRUE)

###Disparity among groups

library("dispRity")

#Assigning species to groups (Andes and Amazon)

Biogeography<-list("Andes" = c("Palicourea_stenosepala","Palicourea_lineata",
                 "Palicourea_stipularis","Palicourea_flavescens",
                 "Palicourea_padifolia","Palicourea_lehmannii",
                 "Palicourea_thyrsiflora","Palicourea_demissa",
                 "Palicourea_amethystina","Palicourea_standleyana",
                 "Palicourea_seemannii","Palicourea_pyramidalis",
                 "Palicourea_luteonivea","Palicourea_sulphurea",
                 "Palicourea_apicata","Palicourea_loxensis",
                 "Palicourea_angustifolia"
#                 ,"Palicourea_flavifolia","Palicourea_bangii","Palicourea_reticulata"
                 ),
     "Amazon" = c("Palicourea_suerrensis","Palicourea_justiciifolia",
                  "Palicourea_ostreophora","Palicourea_quadrifolia",
                  "Palicourea_corymbifera","Palicourea_winkleri",
                  "Palicourea_dichotoma","Palicourea_gracilenta",
                  "Palicourea_obliquinervia","Palicourea_prunifolia",
                  "Palicourea_callithrix","Palicourea_glabra",
                  "Palicourea_didymocarpos","Palicourea_acuminata",
                  "Palicourea_rhodothamna","Palicourea_andina",
                  "Palicourea_croceoides","Palicourea_crocea",
                  "Palicourea_triphylla","Palicourea_lasiantha",
                  "Palicourea_nitidella","Palicourea_macrobotrys",
                  "Palicourea_guianensis","Palicourea_marcgravii",
                  "Palicourea_grandiflora","Palicourea_rigida",
                  "Palicourea_polycephala","Palicourea_egensis",
                  "Palicourea_deflexa","Palicourea_woronovii"),
     "Central" = c("Palicourea_tetragona","Palicourea_domingensis",
                   "Palicourea_pubescens","Palicourea_elata",
                   "Palicourea_correae","Palicourea_berteroana"),
     "Lowland" = c("Palicourea_timbiquensis","Palicourea_acanthacea",
                   "Palicourea_brachiata","Palicourea_glomerulata"),
     "Atlantic_Forest" = c("Psychotria_suterella","Palicourea_brachypoda",
                           "Palicourea_sessilis","Palicourea_divaricata"))

Biogeography2<-list("Andes" = c("Palicourea_stenosepala","Palicourea_lineata",
                 "Palicourea_stipularis","Palicourea_flavescens",
                 "Palicourea_padifolia","Palicourea_lehmannii",
                 "Palicourea_thyrsiflora","Palicourea_demissa",
                 "Palicourea_amethystina","Palicourea_standleyana",
                 "Palicourea_seemannii","Palicourea_pyramidalis",
                 "Palicourea_luteonivea","Palicourea_sulphurea",
                 "Palicourea_apicata","Palicourea_loxensis",
                 "Palicourea_angustifolia"
#                 ,"Palicourea_flavifolia","Palicourea_bangii","Palicourea_reticulata"
                 ),
     "Amazon" = c("Palicourea_suerrensis","Palicourea_justiciifolia",
                  "Palicourea_ostreophora","Palicourea_quadrifolia",
                  "Palicourea_corymbifera","Palicourea_winkleri",
                  "Palicourea_dichotoma","Palicourea_gracilenta",
                  "Palicourea_obliquinervia","Palicourea_prunifolia",
                  "Palicourea_callithrix","Palicourea_glabra",
                  "Palicourea_didymocarpos","Palicourea_acuminata",
                  "Palicourea_rhodothamna","Palicourea_andina",
                  "Palicourea_croceoides","Palicourea_crocea",
                  "Palicourea_triphylla","Palicourea_lasiantha",
                  "Palicourea_nitidella","Palicourea_macrobotrys",
                  "Palicourea_guianensis","Palicourea_marcgravii",
                  "Palicourea_grandiflora","Palicourea_rigida",
                  "Palicourea_polycephala","Palicourea_egensis",
                  "Palicourea_deflexa","Palicourea_woronovii"),
     "Central" = c("Palicourea_tetragona","Palicourea_domingensis",
                   "Palicourea_pubescens","Palicourea_elata",
                   "Palicourea_correae","Palicourea_berteroana"),
     "Lowland" = c("Palicourea_timbiquensis","Palicourea_acanthacea",
                   "Palicourea_brachiata","Palicourea_glomerulata"),
     "Atlantic_Forest" = c("Psychotria_suterella","Palicourea_brachypoda",
                           "Palicourea_sessilis","Palicourea_divaricata"),
     "out" = c("Palicourea_timbiquensis","Palicourea_acanthacea",
               "Palicourea_brachiata","Palicourea_glomerulata",
               "Palicourea_conephoroides","Palicourea_tinctoria",
               "Palicourea_jelskii","Palicourea_allenii",
               "Palicourea_tomentosa","Psychotria_rosea",
               "Palicourea_hazenii","Palicourea_cyanococca",
               "Psychotria_suterella","Palicourea_brachypoda",
               "Palicourea_sessilis","Palicourea_divaricata"))



#########################################
#Regular PCA for dispRity analysis below
#########################################

pc <- prcomp(median_niche_values[,2:5],
             center = TRUE,
             scale. = TRUE)
attributes(pc)
print(pc)
summary(pc)

var <- get_pca_var(pc)
var$contrib

##Create a dataframe for PCA results
pca_df <- data.frame(
  Species = median_niche_values$species,
  PC1 = pc$x[,1],
  PC2 = pc$x[,2],
  PC3 = pc$x[,3],
  PC4 = pc$x[,4]
)
#write.table(pca_df, "Palicourea_niche_PCA.csv", quote = FALSE, sep = "\t", )
#write.table(pca_df, "Palicourea_niche_PCA_with_elev.csv", quote = FALSE, sep = "\t", )
pca_df_rownames <- data.frame(pca_df[,-1], row.names=pca_df[,1])
pca_df_rownames<-pca_df_rownames[,1:2]
pca_df_rownames_pc2<-pca_df_rownames[,2:3]
matrix_pca_df<-as.matrix(pca_df_rownames)
matrix_pca_df_pc2<-as.matrix(pca_df_rownames_pc2)
##Plot the PCA results
#ggplot(pca_df, aes(x = PC1, y = PC2, label = Species)) +
#  geom_point(size = 3) +
#  geom_text(vjust = -1, hjust = 0.5) +
#  theme_minimal()

#fviz_pca_biplot(pc, repel = TRUE,
#                col.var = "#2E9FDF", # Variables color
#                col.ind = "#696969"  # Individuals color
#)

##Disparity among groups using the average squared pairwise distance metric

#Brootstrapping with rarefaction
subsets<-custom.subsets(data=matrix_pca_df, group = Biogeography2)
boot<-boot.matrix(subsets, bootstraps = 100,
            rarefaction = 4)

#Estimating disparity
disparity_rarefied <- dispRity(boot, dimensions=1,metric = function(x) mean(dist(x)^2))
disparity_rarefied
summary(disparity_rarefied)

#Plotting the results
plot(disparity_rarefied, observed =T)
#plot(disparity, observed = list("pch" = 19, col = "blue", cex = 4))

# Measuring the subset overlap
test.dispRity(disparity_rarefied, test = bhatt.coeff, correction = "bonferroni") #in supplements

# Differences between the subsets bootstrapped
test.dispRity(disparity_rarefied, test = wilcox.test, comparisons = "pairwise", correction = "bonferroni",
              conc.quantiles = c(mean, c(95, 5))) #in supplements


#Disparity among groups using the average squared pairwise distance metric and No bootstrap nor rarification
#disparity_data <- dispRity.per.group(matrix_median_niche_values,
#                                     group = Andes_Amazon,
#                                     metric = function(x) mean(dist(x)^2))

#disparity_data
#summary(disparity_data)
#Plotting the results
#plot(disparity_data)

#Testing for a difference between groups
#test.dispRity(disparity_data, test = bhatt.coeff)

## Calculate the disparity of the dataset using dtt.dispRity
dispRity_dtt <- dtt.dispRity(data = matrix_median_niche_values[,2:5], metric = function(x) mean(dist(x)^2),
                             tree = tree, nsim = 100, scale.time=T)

#dispRity_dtt
## Plotting the results
plot(dispRity_dtt)


#########################################################################################
###############################Phylogenetic PCA #########################################
#########################################################################################

##Principal Components Analysis
##Reading tree
astral01 <- read.tree("rev_dendrogram.tre")
tips2delete<-c("Pal_acuminata1","Pal_guianensis2","Pal_quinquepyrena")
tree <- ape::drop.tip(astral01, tips2delete)
plot(tree)

# Update tip labels so that they match the names in the niche data
new_names <- read.csv("pal_name_updates.csv", header = T)
match_indices <- match(tree$tip.label, new_names$old_label)
# Replace names in tree$tip.label vector
tree$tip.label <- ifelse(is.na(match_indices), tree$tip.label, new_names$new_label[match_indices])
plot(tree)

#add Areas to highlight for plotting
areas<-c("L","E","Other","A","E","A","A","Other","C","L","F","Other","E","Other",
         "C","E","E","E","Other","E","A","E","E","F","C","E","C","A",
         "Other","E","L","E","E","E","Other","Other","E","E","A","A","A","A",
         "E","E","E","E","E","A","Other","E","E","C","A","E","Other","A",
         "E","E","A","F","A","A","A","Other","E","A","C","A","L","Other",
         "Other","Other","E","E","E","Other","F") ##clumsy but works. You make it better
median_niche_values$areas<-areas

#convert Species column into row names
median_niche_values_rownames <- data.frame(median_niche_values[,-1], row.names=median_niche_values[,1])

#Convert dataframe to matrix
matrix_median_niche_values<-data.matrix(median_niche_values_rownames)

#Phylogenetic PCA
pPCA<-phyl.pca(tree, matrix_median_niche_values[,1:4])

attributes(pPCA)
print(pPCA)
summary(pPCA)

#Plotting Phylogenetic PCA
#cols<-setNames(palette()[6:1],sort(unique(getStates(tree))))
#par(mar=c(5.1,4.1,0.6,0.6))

obj<-as.princomp(pPCA)

p1<-fviz_screeplot(obj,addlabels=TRUE)
p2<-fviz_pca_var(obj, col.var="contrib",
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE # Avoid text overlapping
)

p3<-fviz_pca_ind(obj,label="none", habillage=median_niche_values$areas,
                 palette= c("deepskyblue","#FFFF00","forestgreen", 
                            "yellowgreen", "darkorange1","grey" ), addEllipses=F, pointsize=3,
                 repel=T, max.overlaps=100)

plotnames <- p3$data$name %>% as.character()
plotnames<-recode(plotnames, "Palicourea_divaricata"="Pal_dicarivata", 
                  "Palicourea_timbiquensis"="Pal_timbiquensis")

mylabels <- sapply(plotnames, function(x) ifelse(is.na(str_locate(x,"Pal_")[1]),
                                                 "", x)) %>% as.vector()
p3<-p3 + geom_text(aes(label = mylabels))
ggpubr::ggarrange(p3)
ggpubr::ggarrange(p2) #in supplements

#phylomorphospace(tree,scores(pali_pca,dim=c(1,2)),
#                 ftype="off",node.by.map=TRUE,bty="n",
#                 node.size=c(0,1.2))
#grid()
#legend("topleft",names(cols),pch=21,pt.bg=cols,horiz=TRUE,
#       bty="n",pt.cex=1.5)

#########################################################################################################
########################### Inflorescence type, climate, and elevation##################################
#########################################################################################################

data<-read.csv("Palicourea_median_climate_data_with_elev_pollination.csv")

pollination<-data$Principal_pollinator
species<-data$species
pollination<- setNames(pollination, species)

elevation<-data$elev
elevation<- setNames(elevation, species)

#bio12<-data$bio12
#bio12<- setNames(bio12, species)
#convert Species column into row names

phylANOVA(tree, pollination, elevation)

p <- ggplot(data, aes(x=Principal_pollinator, y=elev, color=Principal_pollinator)) + # fill=name allow to automatically dedicate a color for each group
  geom_boxplot()
p

#phylANOVA(tree, pollination, bio12)
#PLOT THAT HERE TO SHOW THAT THERE IS THAT TENDENCY BUT THEY CAN ALSO GO DOWN

#one.way <- aov(elev ~ Principal_pollinator, data = data)

#########################################################################################################
########################### Calculating Niche Overlap with nicheRover##################################
#########################################################################################################

#explore_data
aggregate(niche_data_complete[2:5], niche_data_complete[1], median)
#Prepare data
new_names <- as.character(unique(niche_data_complete$species))
new_names_sort<-sort(new_names)
palicourea_split<-split(niche_data_complete, niche_data_complete$species)
for (i in 1:length(palicourea_split)) {
  assign(new_names_sort[i], palicourea_split[[i]])
}

#I created a file with combinations matching the order in the phylogeny
combinations<-read.csv("combinations.csv",h=T)
spp1<-as.character(combinations$spp1)
spp2<-as.character(combinations$spp2)

pairs <- c(1:3003)

for (i in seq_along(pairs)){
  assign(paste('pair_', i, sep = ''), rbind(get(spp1[i]),get(spp2[i])))
}

#Generate parameter draws from the "default" posterior for each species
nsamples <- 1e3
pair_names <- paste0("pair_", 1:3003)
for (i in 1:length(pair_names)){
  assign(paste('niche_data_complete.par',i, sep = ''), tapply(1:nrow(get(pair_names[i])), get(pair_names[i])$species,
                      function(ii) niw.post(nsamples = nsamples, X = get(pair_names[i])[ii,2:5])))
}

#save this precious info
lapply(ls(pattern="niche_data_complete.par[1-9]+"),function(x) save(list=x, file=paste0(x,".RData")))

#open them
file_names <- paste0("niche_data_complete.par", 1:3003,".Rdata")
for(i in 1:length(file_names)) load(file_names[[i]]) 


#uncomment below for graphs with posterior distributions of each climatic variable in species pairs
#par(mar = c(2, 2, .5, .1)+.1, mfrow = c(1,3))
#clrs <- c("blue", "orange")
#niche.par.plot(niche_data_complete.par2, col = clrs, plot.index = 1)
#niche.par.plot(niche_data_complete.par2, col = clrs, plot.index = 2)
#niche.par.plot(niche_data_complete.par2, col = clrs, plot.index = 1:2)
#legend("topleft", legend = names(niche_data_complete.par), fill = clrs)

#Overlap calculation
niche_data.par_names <- paste0("niche_data_complete.par", 1:3003)
for (i in 1:length(niche_data.par_names)){
  assign(paste('over.stat',i, sep = ''), overlap(get(niche_data.par_names[i]), nreps = nsamples, nprob = 1e3, alpha = .95))
}

#save this precious info
lapply(ls(pattern="over.stat[1-9]+"),function(x) save(list=x, file=paste0(x,".RData")))


#open them
file_names <- paste0("over.stat", 1:3003,".Rdata")
for(i in 1:length(file_names)) load(file_names[[i]]) 

####If making plots of prob. overlap##############################################################################
clrs <- c("blue", "orange")
overlap.plot(over.stat1469, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")#for (i in 1:length(over.stat_names)){
#  overlap.plot((get(over.stat_names[2])), col = clrs, mean.cred.col = "turquoise",equal.axis = TRUE,
#              xlab = "Overlap Probability (%) -- Niche Region Size: 95%")
#}
#######################################################################################################

#Clumsy way of removing same species comparison (diagonal)
diagonal<-c(1,78,154,229,303,376,448,519,589,658,726,793,859,924,988,
            1051,1113,1174,1234,1293,1351,1408,1464,1519,1573,1626,1678,1729,1779,
            1828,1876,1923,1969,2014,2058,2101,2143,2184,2224,2263,2301,2338,2374,
            2409,2443,2476,2508,2539,2569,2598,2626,2653,2679,2704,2728,
            2751,2773,2794,2814,2833,2851,2868,2884,2899,2913,2926,2938,
            2949,2959,2968,2976,2983,2989,2994,2998,3001,3003)
over.stat_names <- paste0("over.stat", 1:3003)
over.stat_names_diagonal<-paste0("over.stat",diagonal)
over.stat_names_final<-setdiff(over.stat_names,over.stat_names_diagonal)

##Calculating mean over all iterations for niche overlap
for (i in 1:length(over.stat_names_final)){
  assign(paste('over.mean',i, sep = ''), apply(get(over.stat_names_final[i]), c(1,2), mean)*100)
}

####################################################################################################
##########PLOTTING##################################################################################
####################################################################################################

#Creating dataframe for plotting
all_species4df<-data.frame(spp1,spp2)
#another clumsy solution for removing same species comparison (diagonal)
species4df<-all_species4df[-c(1,78,154,229,303,376,448,519,589,658,726,793,859,924,988,
                              1051,1113,1174,1234,1293,1351,1408,1464,1519,1573,1626,1678,1729,1779,
                              1828,1876,1923,1969,2014,2058,2101,2143,2184,2224,2263,2301,2338,2374,
                              2409,2443,2476,2508,2539,2569,2598,2626,2653,2679,2704,2728,
                              2751,2773,2794,2814,2833,2851,2868,2884,2899,2913,2926,2938,
                              2949,2959,2968,2976,2983,2989,2994,2998,3001,3003),]

#Creating vectors with overlap mean calculated values for each species pair
over.mean_names <- paste0("over.mean", 1:2926)
overlap_mean_sp1to2=c()
for (i in 1:length(over.mean_names)){
  overlap_mean_sp1to2=c(overlap_mean_sp1to2, get(over.mean_names[i]) [2])
}

overlap_mean_sp2to1=c()
for (i in 1:length(over.mean_names)){
  overlap_mean_sp2to1=c(overlap_mean_sp2to1,  get(over.mean_names[i]) [3])
}

overlaps<-data.frame(overlap_mean_sp1to2, overlap_mean_sp2to1)
plotting_df<-cbind(species4df,overlaps)

write.csv(plotting_df, file = "~/Desktop/Palicourea/Manuscript/niche/plotting_df_revised.csv")

##Creating dendrogram from phylogenetic tree for Heatmap
setwd ("~/Desktop/Palicourea/Manuscript/niche")
astral<-"rev_dendrogram.tre"
astral01<-ape::read.tree(astral)
astral01
#tips2delete<-read.table("tips2delete4dendrogram.txt",header=F)
tips2delete<-c("Pal_acuminata1","Pal_guianensis2","Pal_quinquepyrena")
tree <- ape::drop.tip(astral01, tips2delete)
#check label no. is the same as niche data
tree$tip.label
setwd ("~/Desktop/Palicourea/Manuscript/niche/")
write.tree(tree,file="dendrogram.tre")
#Manually check labels on dendrogram to labels on niche data then import
dendro<-ape::read.tree("dendrogram_edited.tre")

####Heatmap######
#set dendrogram
dendrogram<-as.dendrogram(dendro)
#added values for same species to 100 overlap manually
plotting_df<-read.csv("heatmap.csv",header=T)
heatmap<-plotting_df[,2:78]
heatmap<-as.matrix(heatmap)
rownames(heatmap) = c("Palicourea_suerrensis","Palicourea_justiciifolia","Palicourea_ostreophora","Palicourea_quadrifolia",
                      "Palicourea_corymbifera","Palicourea_winkleri","Palicourea_dichotoma","Palicourea_gracilenta",
                      "Palicourea_obliquinervia","Palicourea_prunifolia","Palicourea_callithrix","Palicourea_glabra",
                      "Palicourea_acuminata","Palicourea_andina","Palicourea_didymocarpos","Palicourea_rhodothamna",
                      "Palicourea_triphylla","Palicourea_croceoides","Palicourea_crocea","Palicourea_lasiantha",
                      "Palicourea_nitidella","Palicourea_macrobotrys","Palicourea_guianensis","Palicourea_marcgravii",
                      "Palicourea_grandiflora","Palicourea_rigida","Palicourea_egensis","Palicourea_polycephala",
                      "Palicourea_deflexa","Palicourea_woronovii","Palicourea_bangii","Palicourea_flavifolia",
                      "Palicourea_reticulata","Palicourea_stipularis","Palicourea_flavescens","Palicourea_stenosepala",
                      "Palicourea_lineata","Palicourea_padifolia","Palicourea_lehmannii","Palicourea_thyrsiflora",
                      "Palicourea_demissa","Palicourea_amethystina","Palicourea_standleyana","Palicourea_seemannii",
                      "Palicourea_pyramidalis","Palicourea_luteonivea","Palicourea_sulphurea","Palicourea_apicata",
                      "Palicourea_loxensis","Palicourea_angustifolia","Palicourea_tetragona","Palicourea_domingensis",
                      "Palicourea_pubescens","Palicourea_correae","Palicourea_elata","Palicourea_berteroana",
                      "Palicourea_petiolaris","Palicourea_timbiquensis","Palicourea_acanthacea","Palicourea_brachiata",
                      "Palicourea_glomerulata","Palicourea_conephoroides","Palicourea_tinctoria","Palicourea_jelskii",
                      "Palicourea_allenii","Palicourea_cyanococca","Palicourea_hazenii","Palicourea_tomentosa",
                      "Psychotria_rosea","Psychotria_suterella","Palicourea_brachypoda","Palicourea_sessilis",
                      "Palicourea_divaricata","Palicourea_racemosa","Palicourea_subfusca","Palicourea_topoensis","Palicourea_brevicollis")
colnames(heatmap) = c("Palicourea_suerrensis","Palicourea_justiciifolia","Palicourea_ostreophora","Palicourea_quadrifolia",
                      "Palicourea_corymbifera","Palicourea_winkleri","Palicourea_dichotoma","Palicourea_gracilenta",
                      "Palicourea_obliquinervia","Palicourea_prunifolia","Palicourea_callithrix","Palicourea_glabra",
                      "Palicourea_acuminata","Palicourea_andina","Palicourea_didymocarpos","Palicourea_rhodothamna",
                      "Palicourea_triphylla","Palicourea_croceoides","Palicourea_crocea","Palicourea_lasiantha",
                      "Palicourea_nitidella","Palicourea_macrobotrys","Palicourea_guianensis","Palicourea_marcgravii",
                      "Palicourea_grandiflora","Palicourea_rigida","Palicourea_egensis","Palicourea_polycephala",
                      "Palicourea_deflexa","Palicourea_woronovii","Palicourea_bangii","Palicourea_flavifolia",
                      "Palicourea_reticulata","Palicourea_stipularis","Palicourea_flavescens","Palicourea_stenosepala",
                      "Palicourea_lineata","Palicourea_padifolia","Palicourea_lehmannii","Palicourea_thyrsiflora",
                      "Palicourea_demissa","Palicourea_amethystina","Palicourea_standleyana","Palicourea_seemannii",
                      "Palicourea_pyramidalis","Palicourea_luteonivea","Palicourea_sulphurea","Palicourea_apicata",
                      "Palicourea_loxensis","Palicourea_angustifolia","Palicourea_tetragona","Palicourea_domingensis",
                      "Palicourea_pubescens","Palicourea_correae","Palicourea_elata","Palicourea_berteroana",
                      "Palicourea_petiolaris","Palicourea_timbiquensis","Palicourea_acanthacea","Palicourea_brachiata",
                      "Palicourea_glomerulata","Palicourea_conephoroides","Palicourea_tinctoria","Palicourea_jelskii",
                      "Palicourea_allenii","Palicourea_cyanococca","Palicourea_hazenii","Palicourea_tomentosa",
                      "Psychotria_rosea","Psychotria_suterella","Palicourea_brachypoda","Palicourea_sessilis",
                      "Palicourea_divaricata","Palicourea_racemosa","Palicourea_subfusca","Palicourea_topoensis","Palicourea_brevicollis")

#col_fun = colorRamp2(c(0,50,100), c("#4EB265","#F7F056","#FE9929"))
col_fun = colorRamp2(c(0,50,100), c("#F7F056","#FE9929","#DC050C"))

Heatmap(heatmap, cluster_rows=dendrogram, cluster_columns=dendrogram, column_dend_height = unit(2, "cm"),
        row_dend_width = unit(2, "cm"),column_names_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 5), col= colorRampPalette(brewer.pal(9, "OrRd"))(10))
#

             
#commands for each step above (no loops)
#over.stat <- overlap(niche_data_complete.par300, nreps = 1000, nprob = 1e4, alpha = .95
#over.mean <- apply(over.stat, c(1:2,1), mean)*100
#round(over.mean, 2)
#over.cred <- apply(over.stat*100, c(1:2, 4), quantile, prob = c(.025, .975), na.rm = TRUE)
#round(over.cred[,,,1])
#over.stat7 <- overlap(niche_data_complete.par7, nreps = nsamples, nprob = 1e4, alpha = .95)
#overlap.plot(over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
#             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")



########################################################################################################
########################################################################################################
########################SUBSET HEATMAP##################################################################
########################################################################################################

#heatmap_all<-read_csv("plotting_df.csv")

#keep<-c("Palicourea_justiciifolia","Palicourea_ostreophora","Palicourea_quadrifolia",
#        "Palicourea_corymbifera","Palicourea_winkleri","Palicourea_dichotoma","Palicourea_gracilenta",
#        "Palicourea_obliquinervia","Palicourea_prunifolia","Palicourea_callithrix","Palicourea_glabra",
#        "Palicourea_acuminata","Palicourea_didymocarpos","Palicourea_rhodothamna","Palicourea_triphylla",
#        "Palicourea_croceoides","Palicourea_crocea","Palicourea_lasiantha","Palicourea_nitidella",
#        "Palicourea_grandiflora","Palicourea_rigida","Palicourea_egensis","Palicourea_polycephala",
#        "Palicourea_deflexa","Palicourea_woronovii")



#dendro_subset<-ape::read.tree("dendrogram_subset.tre")
#dendro_subset<-as.dendrogram(dendro_subset)
####Heatmap######

#heatmap_subset<-read.csv("heatmap_subset.csv",header=T)
#added values for same secies to 100 overlap manually
#heatmap_subset<-heatmap_subset[,2:26]
#heatmap_subset<-as.matrix(heatmap_subset)
#rownames(heatmap_subset) = c("Palicourea_justiciifolia","Palicourea_ostreophora","Palicourea_quadrifolia",
#                      "Palicourea_corymbifera","Palicourea_winkleri","Palicourea_dichotoma","Palicourea_gracilenta",
#                      "Palicourea_obliquinervia","Palicourea_prunifolia","Palicourea_callithrix","Palicourea_glabra",
#                      "Palicourea_acuminata","Palicourea_didymocarpos","Palicourea_rhodothamna",
#                      "Palicourea_triphylla","Palicourea_croceoides","Palicourea_crocea","Palicourea_lasiantha",
#                      "Palicourea_nitidella","Palicourea_grandiflora","Palicourea_rigida","Palicourea_egensis","Palicourea_polycephala",
#                      "Palicourea_deflexa","Palicourea_woronovii")
#colnames(heatmap_subset) = c("Palicourea_justiciifolia","Palicourea_ostreophora","Palicourea_quadrifolia",
#                      "Palicourea_corymbifera","Palicourea_winkleri","Palicourea_dichotoma","Palicourea_gracilenta",
#                      "Palicourea_obliquinervia","Palicourea_prunifolia","Palicourea_callithrix","Palicourea_glabra",
#                      "Palicourea_acuminata","Palicourea_didymocarpos","Palicourea_rhodothamna",
#                      "Palicourea_triphylla","Palicourea_croceoides","Palicourea_crocea","Palicourea_lasiantha",
#                      "Palicourea_nitidella","Palicourea_grandiflora","Palicourea_rigida","Palicourea_egensis","Palicourea_polycephala",
#                      "Palicourea_deflexa","Palicourea_woronovii")

#col_fun = colorRamp2(c(0,50,100), c("#4EB265","#F7F056","#FE9929"))
#col_fun = colorRamp2(c(0,50,100), c("#F7F056","#FE9929","#DC050C"))

#Heatmap(heatmap_subset, cluster_rows=dendro_subset, cluster_columns=dendro_subset, column_dend_height = unit(2, "cm"),
#        row_dend_width = unit(2, "cm"),column_names_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 5), col= colorRampPalette(brewer.pal(9, "OrRd"))(10))
#

#heatmap_subset2<-read.csv("heatmap_subsetF.csv",header=T)
#added values for same secies to 100 overlap manually
#heatmap_subset2<-heatmap_subset2[,2:6]
#heatmap_subset2<-as.matrix(heatmap_subset2)
#rownames(heatmap_subset2) = c("Psychotria_suterella","Palicourea_brachypoda","Palicourea_sessilis",
#                              "Palicourea_divaricata","Palicourea_brevicollis")
#colnames(heatmap_subset2) = c("Psychotria_suterella","Palicourea_brachypoda","Palicourea_sessilis",
#                              "Palicourea_divaricata","Palicourea_brevicollis")
#dendro_subsetF<-ape::read.tree("dendrogram_subsetF.tre")
#dendro_subsetF<-as.dendrogram(dendro_subsetF)

#Heatmap(heatmap_subset2, cluster_rows=dendro_subsetF, cluster_columns=dendro_subsetF, column_dend_height = unit(2, "cm"),
#        row_dend_width = unit(2, "cm"),column_names_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 5), col= colorRampPalette(brewer.pal(9, "OrRd"))(10))

#################################################################
##Proportion of taxa by inflorescence type per evaluated clades##
#################################################################

#inflorescence_types<-read.csv("~/Desktop/Submitted/palicourea/Manuscript/area_inflorescence.csv",header=T)
#io <- inflorescence_types[, percentage := V2 / sum(V2) * 100, by=V3]
#                                                                   labels = c("Short", "Medium", "Long", "Extra Long")))

# Use the table() function to get the counts of each species
#Amazon<-filter(inflorescence_types,area=='Amazon')
#Andes<-filter(inflorescence_types,area=='Andes')
#Central<-filter(inflorescence_types,area=='Central') not enough data
#Lowland<-filter(inflorescence_types,area=='Lowland') not enough data
#out<-filter(inflorescence_types,area=='out')

#i0group_counts<-table(inflorescence_types)
#<- inflorescence_types %>%
#  group_by(V3) %>%
#  summarise(percentage = n() / nrow(inflorescence_types))
