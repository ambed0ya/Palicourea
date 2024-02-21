#########################################################################################################
### 1.Extracting and cleaning occurrence and environmental data for large datasets with CANDI pipeline###
# Modified from Tribble et al. 2023.
#########################################################################################################
setwd("~/Desktop/Palicourea/Manuscript/niche")

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
library("phytools")
library("phangorn")
library("phylogram")
library("circlize")
library('RColorBrewer')

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

load("~/Desktop/Palicourea/Manuscript/niche/occurrences_filtered.RData")

coords<-c("longitude","latitude")
coordinates = lapply(occurrences, "[", , coords)
list2env(coordinates,envir=.GlobalEnv)

#download 19 bioclim variables from WorldClim and process into a RasterStack
climate_stack <- get_world_clim()

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

niche_data<-df[4:8]
niche_data$species<-gsub(" ", "_", niche_data$species)
niche_data_complete<-na.omit(niche_data)



#########################################################################################################
########################### 2.Calculating Niche Overlap with nicheRover##################################
#########################################################################################################

#explore_data
aggregate(niche_data_complete[2:5], niche_data_complete[1], mean)
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
overlap.plot(over.stat70, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")

#for (i in 1:length(over.stat_names)){
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
setwd ("~/Desktop/Palicourea/RevBayes/time_stratified/")
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
#Manually check labels on dendrogram to for labels on niche data then import
dendro<-ape::read.tree("dendrogram_edited.tre")

####Heatmap######
#set dendrogram
dendrogram<-as.dendrogram(dendro)
#added values for same secies to 100 overlap manually
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
#over.stat <- overlap(niche_data_complete.par3850, nreps = nsamples, nprob = 1e3, alpha = c(.95,0.99))
#over.mean <- apply(over.stat, c(1:2,4), mean)*100
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


##TO DO: add subsetting commands here

dendro_subset<-ape::read.tree("dendrogram_subset.tre")
dendro_subset<-as.dendrogram(dendro_subset)
####Heatmap######

heatmap_subset<-read.csv("heatmap_subset.csv",header=T)
#added values for same secies to 100 overlap manually
heatmap_subset<-heatmap_subset[,2:26]
heatmap_subset<-as.matrix(heatmap_subset)
rownames(heatmap_subset) = c("Palicourea_justiciifolia","Palicourea_ostreophora","Palicourea_quadrifolia",
                      "Palicourea_corymbifera","Palicourea_winkleri","Palicourea_dichotoma","Palicourea_gracilenta",
                      "Palicourea_obliquinervia","Palicourea_prunifolia","Palicourea_callithrix","Palicourea_glabra",
                      "Palicourea_acuminata","Palicourea_didymocarpos","Palicourea_rhodothamna",
                      "Palicourea_triphylla","Palicourea_croceoides","Palicourea_crocea","Palicourea_lasiantha",
                      "Palicourea_nitidella","Palicourea_grandiflora","Palicourea_rigida","Palicourea_egensis","Palicourea_polycephala",
                      "Palicourea_deflexa","Palicourea_woronovii")
colnames(heatmap_subset) = c("Palicourea_justiciifolia","Palicourea_ostreophora","Palicourea_quadrifolia",
                      "Palicourea_corymbifera","Palicourea_winkleri","Palicourea_dichotoma","Palicourea_gracilenta",
                      "Palicourea_obliquinervia","Palicourea_prunifolia","Palicourea_callithrix","Palicourea_glabra",
                      "Palicourea_acuminata","Palicourea_didymocarpos","Palicourea_rhodothamna",
                      "Palicourea_triphylla","Palicourea_croceoides","Palicourea_crocea","Palicourea_lasiantha",
                      "Palicourea_nitidella","Palicourea_grandiflora","Palicourea_rigida","Palicourea_egensis","Palicourea_polycephala",
                      "Palicourea_deflexa","Palicourea_woronovii")

#col_fun = colorRamp2(c(0,50,100), c("#4EB265","#F7F056","#FE9929"))
col_fun = colorRamp2(c(0,50,100), c("#F7F056","#FE9929","#DC050C"))

Heatmap(heatmap_subset, cluster_rows=dendro_subset, cluster_columns=dendro_subset, column_dend_height = unit(2, "cm"),
        row_dend_width = unit(2, "cm"),column_names_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 5), col= colorRampPalette(brewer.pal(9, "OrRd"))(10))
#

heatmap_subset2<-read.csv("heatmap_subsetF.csv",header=T)
#added values for same secies to 100 overlap manually
heatmap_subset2<-heatmap_subset2[,2:6]
heatmap_subset2<-as.matrix(heatmap_subset2)
rownames(heatmap_subset2) = c("Psychotria_suterella","Palicourea_brachypoda","Palicourea_sessilis",
                              "Palicourea_divaricata","Palicourea_brevicollis")
colnames(heatmap_subset2) = c("Psychotria_suterella","Palicourea_brachypoda","Palicourea_sessilis",
                              "Palicourea_divaricata","Palicourea_brevicollis")
dendro_subsetF<-ape::read.tree("dendrogram_subsetF.tre")
dendro_subsetF<-as.dendrogram(dendro_subsetF)

Heatmap(heatmap_subset2, cluster_rows=dendro_subsetF, cluster_columns=dendro_subsetF, column_dend_height = unit(2, "cm"),
        row_dend_width = unit(2, "cm"),column_names_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 5), col= colorRampPalette(brewer.pal(9, "OrRd"))(10))

#########################################################################################################
################################EXTRA THINGS TO DECIDE IF DELETE OR LEAVE################################
#########################################################################################################



#######################################################
#######################PCA#############################
#######################################################

#pc<-prcomp(pca_data_complete[,-1])
#write.csv(pca_data,file="~/Desktop/test.csv")
#summary(pc)
#autoplot(pc)
#autoplot(pc, data = pca_data_complete,
#         loadings = TRUE, loadings.colour = 'blue',
#         loadings.label = TRUE, loadings.label.size = 3)
#autoplot(pc, data=pca_data_complete, colour="species",legend="none")

final_list <- split( niche_data_complete , f = niche_data_complete$species )

#nrows=lapply(occurrences,function(x) nrow(x))

mean1<-lapply(final_list, function(x) mean(x$bio1))
df1<-list.rbind(mean1)
#mean2<-lapply(final_list, function(x) mean(x$bio2))
#df2<-list.rbind(mean2)  
#mean3<-lapply(final_list, function(x) mean(x$bio3))
#df3<-list.rbind(mean3)
#mean4<-lapply(final_list, function(x) mean(x$bio4))
#df4<-list.rbind(mean4)
mean5<-lapply(final_list, function(x) mean(x$bio5))
df5<-list.rbind(mean5)  
mean12<-lapply(final_list, function(x) mean(x$bio12))
df12<-list.rbind(mean6)
#mean7<-lapply(final_list, function(x) mean(x$bio7))
#df7<-list.rbind(mean7)
#mean8<-lapply(final_list, function(x) mean(x$bio8))
#df8<-list.rbind(mean8)  
#mean9<-lapply(final_list, function(x) mean(x$bio9))
#df9<-list.rbind(mean9)
#mean10<-lapply(final_list, function(x) mean(x$bio10))
#df10<-list.rbind(mean10)
#mean11<-lapply(final_list, function(x) mean(x$bio11))
#df11<-list.rbind(mean11)  
mean12<-lapply(final_list, function(x) mean(x$bio12))
df12<-list.rbind(mean12)
#mean13<-lapply(final_list, function(x) mean(x$bio13))
#df13<-list.rbind(mean13)
#mean14<-lapply(final_list, function(x) mean(x$bio14))
#df14<-list.rbind(mean14)  
#mean15<-lapply(final_list, function(x) mean(x$bio15))
#df15<-list.rbind(mean15)
#mean16<-lapply(final_list, function(x) mean(x$bio16))
#df16<-list.rbind(mean16)
#mean17<-lapply(final_list, function(x) mean(x$bio17))
#df17<-list.rbind(mean17)
#mean18<-lapply(final_list, function(x) mean(x$bio18))
#df18<-list.rbind(mean18)  
#mean19<-lapply(final_list, function(x) mean(x$bio19))
#df19<-list.rbind(mean19) 
#dfs_total<-cbind(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,df17,df18,df19)
dfs_total<-cbind(df1,df5,df6,df12)
dfs<-as.data.frame(dfs_total)

my_species<-read.table("~/Desktop/Palicourea/species.txt",h=T)
dfs_final<-cbind(my_species,dfs)
#colnames(dfs_final)<-c("species","bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11",
#                 "bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")
colnames(dfs_final)<-c("species","bio1","bio5","bio6",
                       "bio12")
#list2env(final_list,envir=.GlobalEnv)
#use lapply to do PCA on each data fram of list
#pca<-lapply(final_list, function(x) prcomp(x[,-1]))
#list2env(pca,envir=.GlobalEnv)
#autoplot(Palicourea_acanthacea)
PCA_ALL<-prcomp(dfs_final[,-1],TRUE)
#PCA_PREC<-prcomp(dfs_final[,2:5],TRUE)
#PCA_TEMP<-prcomp(dfs_final[,],TRUE)
#autoplot(PCA_PREC)
#autoplot(PCA_TEMP)
autoplot(PCA_ALL)

#######################################################
#######Niche stochastic character mapping##############
#######################################################
#read tree and PCA values
tree <- read.tree("rev_starting.tre")
niche_pca<-PCA_ALL$x[,1]

#Manually remove unwanted values and NA for species that do not have values. I also need to change tip labels. Although this can be done when conducting species tree inference
tips<- tree$tip.label


#do not use#fit<-fastAnc(tree,niche_pca,vars=T,CI=T) or contmap. Shee which one

######################
## MAXENT MODELLING ##
######################

#the jar file evaluates 
#the performance of models built with different subsets of
#the environmental variables so each variable can be given a contribution score
#jar file doesn't come with maxent so it has to be dowloaded separately
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')

results <- model_niches(occ_data = occ_data1, clim_data = climate_stack, reps = 2, reptype = bootstrap, output_dir = "~/Desktop/Palicourea/maxent_output/")
save(results, file = "~/Desktop/Palicourea/maxent_output/maxent_results14.RData")
load("~/Desktop/Palicourea/maxent_output/maxent_results14.RData")

#### Calculate the 'optimal' climate values for each bioclim variable ####

max_wc <- as.data.frame(matrix(nrow = length(results[[2]]),
                               ncol = length(names(final_stack))))
colnames(max_wc) <- names(final_stack)
rownames(max_wc) <- names(results[[2]])

for (i in 1:nrow(max_wc)) {
  print(paste0("working on ", rownames(max_wc)[i], ", species number ", i))
  coords_of_max <- data.frame(long = rasterToPoints(results[[2]][[i]])[which.max( rasterToPoints(results[[2]][[i]])[,3]),c(1,2)][1],
                              lat = rasterToPoints(results[[2]][[i]])[which.max( rasterToPoints(results[[2]][[i]])[,3]),c(1,2)][2])
  max_wc[i,] <- raster::extract(x = final_stack, y = SpatialPoints(coords_of_max))
  
}

#write.csv(max_wc, file = "~/Desktop/Palicourea/maxent_output/max_wc1.csv")
#max_wc <- read.csv("~/Desktop/Palicourea/maxent_output/max_wc1.csv", row.names = 1)
write.table(max_wc, file = "~/Desktop/Palicourea/maxent_output/maxent_output_max_wc1.csv", append = TRUE, sep=",",row.names = T, col.names = T)

#### Compare average world clim values pre modeling to the optimum values from maxent ####

mean_premodeling <- as.data.frame(matrix(nrow = length(occ_data1),
                                         ncol = length(names(final_stack))))
rownames(mean_premodeling) <- names(occ_data1)
colnames(mean_premodeling) <- names(final_stack)

median_premodeling <- as.data.frame(matrix(nrow = length(occ_data1),
                                         ncol = length(names(final_stack))))
rownames(median_premodeling) <- names(occ_data1)
colnames(median_premodeling) <- names(final_stack)

for (i in 1:nrow(mean_premodeling)){
  print(paste0("working on ", rownames(mean_premodeling)[i], ", species number ", i))
  coords <- data.frame(long = occ_data1[[i]]$longitude, 
                       lat = occ_data1[[i]]$latitude)
  clim_values <- raster::extract(x = final_stack, y = SpatialPoints(coords))

  mean_premodeling[i, ] <- (clim_values)
  median_premodeling[i,] <- matrixStats::colMedians(clim_values)
}

#### Boxplots of differences ####
load("~/Desktop/Palicourea/maxent_output/maxent_results1.RData")
# mean occ vs maxent
combined <- cbind(max_wc, mean_premodeling)
colnames(combined) <- c(paste0(colnames(max_wc), "_max_wc"),
                        paste0(colnames(mean_premodeling), "_mean_premodeling"))
combined_matrix <- log(as.matrix(combined) + 422)
combined_mod <- as.data.frame(combined_matrix)


pdf("~/Desktop/maxent_vs_mean_hist1.pdf", height = 10, width = 15)
boxplot(combined_mod[,order(colnames(combined_mod))],
        yaxt = "n", xaxt = "n", axes = T)
axis(1, las = 3, cex.axis=0.75, 
     at = seq(from = 1.5, to = 37.5, by = 2), 
     labels = colnames(max_wc))
rect(xleft = seq(from = 0.5, to = 37.5, by = 2),
     xright = seq(from = 2.5, to = 38.5, by = 2),
     ybottom = rep(4, times = 19),
     ytop = rep(10, times = 19), 
     col = c(rep(c(NA, "grey88"), times = 9), NA),
     lty = 0)
boxplot(combined_mod[,order(colnames(combined_mod))],
        col = rep(c("#EF8A62", "#67A9CF"), times = 19),
        pch = 19, yaxt = "n", xaxt = "n", add = TRUE)
dev.off()

# median occ vs maxent
combined <- cbind(max_wc, median_premodeling)
colnames(combined) <- c(paste0(colnames(max_wc), "_max_wc"),
                        paste0(colnames(median_premodeling), "_median_premodeling"))
combined_matrix <- log(as.matrix(combined) + 422)
combined_mod <- as.data.frame(combined_matrix)


pdf("~/Desktop/Palicourea/maxent_output/maxent_vs_median_hist1.pdf", height = 10, width = 15)
boxplot(combined_mod[,order(colnames(combined_mod))],
        yaxt = "n", xaxt = "n", axes = T)
axis(1, las = 3, cex.axis=0.75, 
     at = seq(from = 1.5, to = 37.5, by = 2), 
     labels = colnames(max_wc))
rect(xleft = seq(from = 0.5, to = 37.5, by = 2),
     xright = seq(from = 2.5, to = 38.5, by = 2),
     ybottom = rep(4, times = 19),
     ytop = rep(10, times = 19), 
     col = c(rep(c(NA, "grey88"), times = 9), NA),
     lty = 0)
boxplot(combined_mod[,order(colnames(combined_mod))],
        col = rep(c("#EF8A62", "#67A9CF"), times = 19),
        pch = 19, yaxt = "n", xaxt = "n", add = TRUE)
dev.off()

#### scatterplots of differences ####

# mean occ vs maxent

vars <- colnames(max_wc)
vars <- vars[c(1,12:19,2:11)]
sample_size <- unlist(lapply(occ_data1, nrow))
rbPal <- colorRampPalette(c('red','blue'))
color <- rbPal(10)[as.numeric(cut(log(sample_size), breaks = 10))]

pdf("~/Desktop/Palicourea/maxent_output/maxent_vs_mean_scatterplots_col1.pdf")
par(mfrow = c(2,2))
for (i in 1:length(vars)) {
  min <- min(c(max_wc[,vars[i]],
               mean_premodeling[,vars[i]]), na.rm = T)
  max <- max(c(max_wc[,vars[i]],
               mean_premodeling[,vars[i]]), na.rm = T)
  diff <- max - min
  if (diff < 200) ticby <- 10
  if (diff > 200 & diff < 1000) ticby <- 50
  if (diff > 1000 & diff < 5000) ticby <- 200
  if (diff > 5000) ticby <- 500
  plot(max_wc[,vars[i]], mean_premodeling[,vars[i]], 
       ylim=c(min,max), xlim = c(min, max), pch = 19,
       xlab = "MaxEnt estimated optimum", ylab = "Mean from occurrence points",
       col = alpha(color, 0.6), axes = F, main = vars[i])
  axis(side = 1, at = seq(from = round(min, digits = -1) , 
                          to = round(max, digits = -1), by = ticby))
  axis(side = 2, at = seq(from = round(min, digits = -1) , 
                          to = round(max, digits = -1), by = ticby))
  abline(a = 0, b = 1, col = "black", lty = 2, lwd = 2)
}
dev.off()

# median occ vs maxent

vars <- colnames(max_wc)
vars <- vars[c(1,12:19,2:11)]

pdf("~/Desktop/maxent_vs_median_scatterplots.pdf")
par(mfrow = c(2,2))
for (i in 1:length(vars)) {
  min <- min(c(max_wc[,vars[i]],
               median_premodeling[,vars[i]]), na.rm = T)
  max <- max(c(max_wc[,vars[i]],
               median_premodeling[,vars[i]]), na.rm = T)
  diff <- max - min
  if (diff < 200) ticby <- 10
  if (diff > 200 & diff < 1000) ticby <- 50
  if (diff > 1000 & diff < 5000) ticby <- 200
  if (diff > 5000) ticby <- 500
  plot(max_wc[,vars[i]], median_premodeling[,vars[i]], 
       ylim=c(min,max), xlim = c(min, max), pch = 19,
       xlab = "MaxEnt estimated optimum", ylab = "median from occurrence points",
       col = alpha("black", 0.6), axes = F, main = vars[i])
  axis(side = 1, at = seq(from = round(min, digits = -1) , 
                          to = round(max, digits = -1), by = ticby))
  axis(side = 2, at = seq(from = round(min, digits = -1) , 
                          to = round(max, digits = -1), by = ticby))
  abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
}
dev.off()

