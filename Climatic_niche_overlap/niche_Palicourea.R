#########################################################################################################
### 1.Extracting and cleaning occurrence and environmental data for large datasets with CANDI pipeline###
# Modified from Tribble et al. 2023.
#########################################################################################################
#setwd("~/Desktop/Submitted/Palicourea/Manuscript/niche")

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
library("sp")
library("geodata")

#source CANDI functions
source("~/Repos/Palicourea/Climatic_niche_overlap/modified_candi/get_occ_records.R")
source("~/Repos/Palicourea/Climatic_niche_overlap/modified_candi/get_both_occ_records.R")
source("~/Repos/Palicourea/Climatic_niche_overlap/modified_candi/get_bien_occ_records.R")
source("~/Repos/Palicourea/Climatic_niche_overlap/modified_candi/get_gbif_occ_records.R")
source("~/Repos/Palicourea/Climatic_niche_overlap/modified_candi/remove_dup_locs.R")
source("~/Repos/Palicourea/Climatic_niche_overlap/modified_candi/remove_ocean_points.R")
source("~/Repos/Palicourea/Climatic_niche_overlap/modified_candi/remove_perf_0_90_180.R")
source("~/Repos/Palicourea/Climatic_niche_overlap/modified_candi/remove_lessthan.R")
source("~/Repos/Palicourea/Climatic_niche_overlap/modified_candi/remove_null_items.R")
source("~/Repos/Palicourea/Climatic_niche_overlap/modified_candi/get_world_clim.R")
source("~/Repos/Palicourea/Climatic_niche_overlap/modified_candi/make_corr_matrix.R")
source("~/Repos/Palicourea/Climatic_niche_overlap/modified_candi/remove_corr_variables.R")
source("~/Repos/Palicourea/Climatic_niche_overlap/modified_candi/remove_points_outside_nat_range.R")
source("~/Repos/Palicourea/Climatic_niche_overlap/modified_candi/model_niches.R")
source("~/Repos/Palicourea/Climatic_niche_overlap/modified_candi/get_world_clim.R")
source("~/Repos/Palicourea/Climatic_niche_overlap/modified_candi/trim_to_shapefile.R")

## OCCURRENCE DATA PREP PATHWAY ##
##################################

#load data# 
#world map
#big_map <- rgdal::readOGR("gadm28_adm0/gadm28_adm0.shp")
big_map <- st_read("~/Bedoya Dropbox/Bedoya_Research_Group/Submitted/palicourea/Manuscript/niche/gadm28_adm0")
big_map_sp <- as(big_map, "Spatial")
#species range data
#read_csv("Palicourea_native_ranges.csv") %>%
#  dplyr::select(species, range) -> all_native_ranges
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
my_species12<-c("Psychotria suterella","Palicourea glabra","Palicourea jelskii","Palicourea nitidella")
my_species13<-c("Palicourea acetosoides","Palicourea aeneofusca","Palicourea aetantha", "Palicourea alagoana", "Palicourea albaniana","Palicourea albiflora","Palicourea albocaerulea")
my_species14<-c("Palicourea amapaensis","Palicourea laxivenulosa","Palicourea amorimii","Palicourea amplissima","Palicourea anceps","Palicourea andaluciana","Palicourea anderssoniana")
my_species15<-c("Palicourea andrei","Palicourea antioquiana","Palicourea antisanana","Palicourea aschersoniana","Palicourea atlantica","Palicourea azulina","Palicourea beachiana")
my_species16<-c("Palicourea bella","Palicourea boqueronensis","Palicourea boraginoides","Palicourea botryocephala","Palicourea breedlovei","Palicourea caerulea","Palicourea calidicola")
my_species17<-c("Palicourea calophylla","Palicourea calothyrsus","Palicourea canarina","Palicourea candida","Palicourea caprifoliacea","Palicourea cauligera","Palicourea cenepensis")
my_species18<-c("Palicourea chignul","Palicourea colorata","Palicourea coriacea","Palicourea diminuta","Palicourea ernestii","Palicourea eurycarpa","Palicourea fanshawei")
my_species19<-c("Palicourea flaviflora","Palicourea frontinoensis","Palicourea fuchsioides","Palicourea fulgens","Palicourea garciae","Palicourea gardenioides","Palicourea gemmiflora")
my_species20<-c("Palicourea graciliflora","Palicourea grandiceps","Palicourea grandifructa","Palicourea grandistipula","Palicourea harlingii","Palicourea heterochroma","Palicourea huampamiensis")
my_species21<-c("Palicourea hypochlorina","Palicourea jambosioides","Palicourea lachnantha","Palicourea lasiorrhachis","Palicourea liogieri","Palicourea lobbii","Palicourea locuples")
my_species22<-c("Palicourea longiflora","Palicourea macarthurorum","Palicourea macbridei","Palicourea mansoana","Palicourea mediocris","Palicourea minutiflora","Palicourea nigricans")
my_species23<-c("Palicourea officinalis","Palicourea palenquensis","Palicourea patens","Palicourea perquadrangularis","Palicourea persearum","Palicourea rigidifolia")
my_species24<-c("Palicourea salicifolia","Palicourea skotakii","Palicourea sodiroi","Palicourea solitudinum","Palicourea sopkinii","Palicourea spathacea","Palicourea stenostachya")
my_species25<-c("Palicourea subspicata","Palicourea sucllii","Palicourea thermydri","Palicourea trianae","Palicourea ulloana","Palicourea umbelliformis","Palicourea valerioana")
my_species26<-c("Palicourea vesiculifera","Palicourea wolffiae","Palicourea yneziae","Psychotria bertieroides","Psychotria bolivarensis","Palicourea bracteocardia")
my_species27<-c("Psychotria bremekampiana","Psychotria everardii","Psychotria hebeclada","Psychotria hemicephaelis","Psychotria limitanea","Psychotria lindenii")
my_species28<-c("Psychotria luxurians","Psychotria marginata","Psychotria muscosa","Psychotria oblonga","Psychotria oinochrophylla","Psychotria panamensis")
my_species29<-c("Psychotria paniculata","Psychotria patentinervia","Psychotria phyllocalymma","Psychotria prancei","Psychotria pseudinundata")
my_species30<-c("Psychotria rhombibractea","Psychotria rhytidocarpa","Psychotria ruelliifolia","Psychotria spathicalyx","Psychotria stachyoides","Psychotria stipulosa")
my_species31<-c("Psychotria ulviformis","Psychotria urceolata","Psychotria variegata","Psychotria venulosa","Carapichea urniformis")


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
occ_data13 <- get_occ_records(species = my_species13, database = "both")
occ_data14 <- get_occ_records(species = my_species14, database = "both")
occ_data15 <- get_occ_records(species = my_species15, database = "both")
occ_data16 <- get_occ_records(species = my_species16, database = "both")
occ_data17 <- get_occ_records(species = my_species17, database = "both")
occ_data18 <- get_occ_records(species = my_species18, database = "both")
occ_data19 <- get_occ_records(species = my_species19, database = "both")
occ_data20 <- get_occ_records(species = my_species20, database = "both")
occ_data21 <- get_occ_records(species = my_species21, database = "both")
occ_data22 <- get_occ_records(species = my_species22, database = "both")
occ_data23 <- get_occ_records(species = my_species23, database = "both")
occ_data24 <- get_occ_records(species = my_species24, database = "both")
occ_data25 <- get_occ_records(species = my_species25, database = "both")
occ_data26 <- get_occ_records(species = my_species26, database = "both")
occ_data27 <- get_occ_records(species = my_species27, database = "both")
occ_data28 <- get_occ_records(species = my_species28, database = "both")
occ_data29 <- get_occ_records(species = my_species29, database = "both")
occ_data30 <- get_occ_records(species = my_species30, database = "both")
occ_data31 <- get_occ_records(species = my_species31, database = "both")

occurrences<-c(occ_data1,occ_data2,occ_data3,occ_data4,occ_data5,occ_data6,occ_data7,
               occ_data8,occ_data9,occ_data10,occ_data11,occ_data12,occ_data13,occ_data14,
               occ_data15,occ_data16,occ_data17,occ_data18,occ_data19,occ_data20,occ_data21,
               occ_data22,occ_data23,occ_data24,occ_data25,occ_data26,occ_data27,occ_data28,
               occ_data29,occ_data30,occ_data31)
save(occurrences, file = "~/Repos/Palicourea/Climatic_niche_overlap/occurrences_raw_allspp.RData")

#data cleaning
occurrences <- remove_dup_locs(occurrences)
occurrences <- remove_ocean_points(occurrences, world_map = big_map_sp)
occurrences <- remove_perf_0_90_180(occurrences)
#occ_data <- remove_points_outside_nat_range(df = occurrences, 
#                                           botan_map = kew_map_level_2, 
#                                           nat_range = all_native_ranges)
occurrences <- remove_lessthan(occurrences, n = 10)
occurrences <- remove_null_items(occurrences)

save(occurrences, file = "~/Repos/Palicourea/Climatic_niche_overlap/occurrences_filtered_allspp.RData")

## CLIMATE DATA PREP PATHWAY ##
###############################

load("~/Repos/Palicourea/Climatic_niche_overlap/occurrences_filtered_allspp.RData")

coords<-c("longitude","latitude")
coordinates = lapply(occurrences, "[", , coords)
list2env(coordinates,envir=.GlobalEnv)

#download 19 bioclim variables from WorldClim and process into a RasterStack
climate_stack <- get_world_clim()

#If the above does not work then download manually
setwd("~/Repos/Palicourea/Climatic_niche_overlap/")
#makes pdf of correlation matrix saved to working directory to pick out correlated variables
matrix <- make_corr_matrix(occurrences = occurrences, environment_data = climate_stack, abs_highlight = 0.8)

#select variables to remove
bad_vars <- c("wc2.1_5m_bio_2", "wc2.1_5m_bio_3", "wc2.1_5m_bio_6", "wc2.1_5m_bio_7",
              "wc2.1_5m_bio_8", "wc2.1_5m_bio_9","wc2.1_5m_bio_10","wc2.1_5m_bio_11",
              "wc2.1_5m_bio_13", "wc2.1_5m_bio_14", "wc2.1_5m_bio_15","wc2.1_5m_bio_16","wc2.1_5m_bio_17",
              "wc2.1_5m_bio_18","wc2.1_5m_bio_19") #test

#remove variables selected by user
uncorr_stack <- remove_corr_variables(environment_data = climate_stack, variables_to_be_removed = bad_vars) 
final_stack <- uncorr_stack
#final_stack <- climate_stack


##EXTRACTING BIOCLIM VALUES FROM OCCURRENCE POINTS##
####################################################
#points <- lapply(coordinates, function(x) SpatialPoints(x, proj4string = final_stack@crs))
points <- lapply(coordinates, function(x) {
  sp::SpatialPoints(x, proj4string = sp::CRS(terra::crs(final_stack)))
})
points_v <- lapply(points, terra::vect)
values <- lapply(points_v, function(p) terra::extract(final_stack, p))
#values <- lapply(points, function(x) extract(final_stack,x))

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

niche_data<-df[3:9]

niche_data$species<-gsub(" ", "_", niche_data$species)
niche_data_complete<-na.omit(niche_data)
write.csv(niche_data_complete, "Palicourea_raw_climate_data_with_elev_final.csv", row.names = F)

##Calculation of median value for climatic variables for each species
median_niche_values<-aggregate(niche_data_complete[3:7], niche_data_complete[1], median)

write.csv(median_niche_values, "Palicourea_median_climate_data_with_elev_final.csv", row.names=F)
#median_niche_values$species<-as.factor(median_niche_values$species)
#load data

median_niche_values<-read.csv("Palicourea_median_climate_data_clades_final.csv")
#########################################################################################
##################DISPARITY ACROSS GROUPS AND THROUGH TIME###############################
#########################################################################################
#convert Species column into row names
median_niche_values_rownames <- data.frame(median_niche_values[,-1], row.names=median_niche_values[,1])

#Convert dataframe to matrix
matrix_median_niche_values<-data.matrix(median_niche_values_rownames)
#Calculate and plot DTT across variables
#par( mfrow= c(5,1) )


###Disparity among groups

library("dispRity")

#Assigning species to groups


Biogeography2<-list("Andes" = c(#"Palicourea_padifolia",
                   "Palicourea_lehmannii",
                   "Palicourea_thyrsiflora","Palicourea_demissa",
                   "Palicourea_frontinoensis", #"Palicourea_angustifolia",
                   "Palicourea_thermydri","Palicourea_stipularis",
                   "Palicourea_flavescens","Palicourea_stenosepala",
                   "Palicourea_albiflora","Palicourea_macbridei",
                   "Palicourea_lineata","Palicourea_amethystina",
                   "Palicourea_heterochroma","Palicourea_anceps",
                   "Palicourea_perquadrangularis",#"Palicourea_albocaerulea",
                   "Palicourea_antioquiana","Palicourea_acetosoides",
                   #"Palicourea_spathacea","Palicourea_skotakii",
                   #"Palicourea_bella",
                   "Palicourea_trianae",
                   #"Palicourea_seemannii","Palicourea_amplissima",
                   #"Palicourea_grandistipula",
                   "Palicourea_sopkinii",
                   "Palicourea_standleyana", "Palicourea_candida",
                   "Palicourea_andrei","Palicourea_sodiroi",
                   "Palicourea_caprifoliacea","Palicourea_graciliflora",
                   "Palicourea_lasiorrhachis","Palicourea_pyramidalis",
                   #"Palicourea_chignul",
                   "Palicourea_lobbii",
                   "Palicourea_canarina","Palicourea_gemmiflora",
                   "Palicourea_luteonivea","Palicourea_harlingii",
                   "Palicourea_calothyrsus","Palicourea_anderssoniana",
                   "Palicourea_garciae","Palicourea_andaluciana",
                   #"Palicourea_salicifolia",
                   "Palicourea_sulphurea",
                   "Palicourea_apicata","Palicourea_sucllii",
                   "Palicourea_azulina","Palicourea_loxensis",
                   "Palicourea_ulloana","Palicourea_quinquepyrena",
                   "Palicourea_rigidifolia","Palicourea_fuchsioides",
                   "Palicourea_albaniana"
#                 ,"Palicourea_flavifolia","Palicourea_bangii","Palicourea_reticulata"
                 ),
     "Amazon" = c("Palicourea_lupulina",#"Palicourea_suerrensis",
                  "Psychotria_lindenii","Palicourea_ostreophora",
                  "Palicourea_quadrifolia","Palicourea_corymbifera",
                  "Palicourea_coriacea","Palicourea_winkleri",
                  "Palicourea_dichotoma","Palicourea_gracilenta",
                  "Psychotria_paniculata","Palicourea_amorimii",
                  "Psychotria_bremekampiana","Psychotria_prancei",
                  "Palicourea_prunifolia","Psychotria_rhombibractea",
                  "Psychotria_variegata","Palicourea_glabra",
                  "Psychotria_ulviformis","Palicourea_callithrix",
                  "Palicourea_rhodothamna","Palicourea_huampamiensis",
                  "Palicourea_didymocarpos","Palicourea_acuminata",
                  #"Palicourea_boraginoides","Palicourea_andina",
                  "Palicourea_diminuta","Palicourea_triphylla",
                  "Palicourea_calophylla",#"Palicourea_subspicata",
                  "Palicourea_croceoides","Palicourea_crocea",
                  "Palicourea_macrobotrys","Palicourea_guianensis",
                  "Palicourea_mansoana","Palicourea_lasiantha",
                  "Palicourea_lachnantha","Palicourea_officinalis",
                  "Palicourea_nitidella","Palicourea_nigricans",
                  "Palicourea_macarthurorum","Palicourea_grandiflora",
                  "Palicourea_amapaensis","Palicourea_marcgravii",
                  "Palicourea_aeneofusca","Palicourea_longiflora",
                  "Palicourea_rigida",#"Palicourea_petiolaris",
                  #"Palicourea_aschersoniana",
                  "Psychotria_everardii",
                  "Psychotria_bertieroides",
                  #"Palicourea_laxivenulosa",
                  "Psychotria_venulosa","Palicourea_deflexa",
                  "Palicourea_polycephala","Palicourea_egensis",
                  "Psychotria_limitanea","Palicourea_woronovii",
                  "Psychotria_oinochrophylla"),
     "Central" = c("Palicourea_umbelliformis","Palicourea_calidicola",
                   "Palicourea_eurycarpa","Palicourea_mediocris",
                   "Palicourea_breedlovei","Palicourea_beachiana",
                   "Palicourea_grandifructa","Palicourea_tetragona",
                   "Palicourea_persearum","Psychotria_hebeclada",
                   "Palicourea_pubescens", "Palicourea_elata",
                   "Palicourea_correae","Palicourea_gardenioides",
                   "Palicourea_domingensis","Psychotria_luxurians",
                   "Psychotria_berteroana"),
     "Amazon2" = c("Psychotria_oblonga","Psychotria_muscosa",
                   #"Palicourea_glomerulata",
                   "Psychotria_hemicephaelis",
                   "Palicourea_fanshawei","Palicourea_aetantha",
                   "Carapichea_urniformis","Psychotria_bolivarensis",
                   "Palicourea_yneziae","Palicourea_grandiceps",
                   "Palicourea_hypochlorina","Palicourea_conephoroides",
                   "Palicourea_flaviflora", "Palicourea_cenepensis",
                   "Palicourea_ernestii",#"Palicourea_locuples",
                   #"Palicourea_antisanana","Palicourea_caerulea",
                   #"Palicourea_brachiata","Palicourea_solitudinum",
                   #"Palicourea_timbiquensis","Palicourea_acanthacea",
                   #"Palicourea_palenquensis","Palicourea_allenii",
                   #"Palicourea_tinctoria",
                   #"Palicourea_wolffiae",
                   #"Palicourea_jelskii",
                   #"Palicourea_cauligera",
                   "Palicourea_botryocephala","Palicourea_tomentosa",
                   "Palicourea_bracteocardia","Palicourea_colorata",
                   "Psychotria_urceolata"
                   #"Palicourea_vesiculifera",
                   #"Palicourea_cyanococca","Palicourea_hazenii"
                   ),
     "Atlantic_Forest" = c("Psychotria_stachyoides","Psychotria_ruelliifolia",
                           "Palicourea_sessilis","Psychotria_rhytidocarpa",
                           "Palicourea_fulgens","Psychotria_spathicalyx",
                           "Psychotria_patentinervia","Psychotria_brachypoda",
                           "Psychotria_suterella","Palicourea_jambosioides",
                           "Palicourea_divaricata","Palicourea_atlantica",
                           #"Palicourea_alagoana",
                           "Psychotria_phyllocalymma"),
     "unstructured" = c("Palicourea_stenostachya","Palicourea_racemosa",
               "Palicourea_subfusca","Palicourea_minutiflora",
               "Palicourea_brevicollis","Palicourea_valerioana",
               "Palicourea_topoensis","Palicourea_boqueronensis"))



#########################################
#Regular PCA for dispRity analysis below
#########################################
pc <- prcomp(median_niche_values[,3:6],
             center = TRUE,
             scale. = TRUE)
#pc <- prcomp(median_niche_values[,3:12],
#             center = TRUE,
#             scale. = TRUE)
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
#pca_df_rownames_pc1<-pca_df_rownames[,1:2]
#pca_df_rownames_pc2<-pca_df_rownames[,2:3]
matrix_pca_df<-as.matrix(pca_df_rownames[, "PC1", drop = FALSE])
matrix_pca_df_pc2<-as.matrix(pca_df_rownames[, "PC2", drop = FALSE])
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
subsets<-custom.subsets(data=matrix_pca_df_pc2, group = Biogeography2)
boot<-boot.matrix(subsets, bootstraps = 1000,
            rarefaction = 4)

#Estimating disparity
disparity_rarefied <- dispRity(boot, dimensions=1,metric = function(x) mean(dist(x)^2))
disparity_rarefied
summary(disparity_rarefied)

#Plotting the results
#dev.off()
plot(disparity_rarefied, observed =T)

#plot(disparity, observed = list("pch" = 19, col = "blue", cex = 4))

# Testing for the subset overlap
test.dispRity(disparity_rarefied, test = bhatt.coeff, correction = "bonferroni") #in supplements


#########################################################################################
###############################Phylogenetic PCA #########################################
#########################################################################################

##Reading tree
astral01 <- read.tree("../Resubmission/rev_dendrogram.tre")
#removing taxa for which there isnt enought climate data, or those that have moved outside of the areas where
#the clade originated
tips2delete<-c("Pal_guianensis2","Pal_demissa2","Pal_angustifolia2","Pal_obliquinervia",
               "Pal_wolffiae","Pal_cauligera","Pal_vesiculifera","Pal_alagoana",
               "Pal_laxivenulosa","Pal_padifolia","Pal_angustifolia1",
               "Pal_albocaerulea","Pal_spathacea","Pal_skotakii",
               "Pal_bella","Pal_seemannii","Pal_amplissima",
               "Pal_grandistipula","Pal_chignul","Pal_salicifolia",
               "Pal_flavifolia","Pal_bangii","Pal_reticulata",
               "Pal_suerrensis","Pal_boraginoides","Pal_andina",
               "Pal_subspicata","Pal_petiolaris","Pal_aschersoniana",
               "Pal_glomerulata","Pal_locuples","Pal_antisanana","Pal_caerulea",
               "Pal_brachiata1","Pal_solitudinum","Pal_timbiquensis",
               "Pal_acanthacea","Pal_palenquensis","Pal_allenii",
               "Pal_tinctoria","Pal_wolffiae","Pal_jelskii",
               "Pal_cauligera","Pal_cyanococca","Pal_hazenii",
               "Pal_alagoana", "Psy_stipulosa", "Psy_pseudinundata")
tree <- ape::drop.tip(astral01, tips2delete)
plot(tree)


# Update tip labels so that they match the names in the niche data
new_names <- read.csv("~/Repos/Palicourea/Climatic_niche_overlap/pal_name_updates.csv", header = T)
match_indices <- match(tree$tip.label, new_names$old_label)
# Replace names in tree$tip.label vector
tree$tip.label <- ifelse(is.na(match_indices), tree$tip.label, new_names$new_label[match_indices])
plot(tree)

write.tree(tree, file = "~/Desktop/tree.tre")
#median_niche_values$areas<-areas
#Check that tree and data have same species
tree_species <- tree$tip.label
data_species <- unique(median_niche_values$species)               


data_not_in_tree <- setdiff(data_species, tree_species)
length(data_not_in_tree)
head(data_not_in_tree)

tree_not_in_data <- setdiff(tree_species, data_species)
length(tree_not_in_data)
head(tree_not_in_data)

#convert Species column into row names
median_niche_values_rownames <- data.frame(median_niche_values[,-1], row.names=median_niche_values[,1])

#Convert dataframe to matrix
matrix_median_niche_values<-data.matrix(median_niche_values_rownames)

#Phylogenetic PCA
pPCA<-phyl.pca(tree, matrix_median_niche_values[,2:5])

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

p3<-fviz_pca_ind(obj,label="none", habillage=median_niche_values$area,
                 palette= c("deepskyblue","#FFFF00","forestgreen", "purple",
                            "yellowgreen","grey","black" ), addEllipses=F, pointsize=3,
                 repel=T, max.overlaps=100)


plotnames <- p3$data$name %>% as.character()
plotnames<-recode(plotnames, "Palicourea_divaricata"="Pal_dicarivata", 
                  "Palicourea_timbiquensis"="Pal_timbiquensis", "Palicourea_acanthacea"="Pal_acanthacea")

mylabels <- sapply(plotnames, function(x) ifelse(is.na(str_locate(x,"Pal_")[1]),
                                                 "", x)) %>% as.vector()
p3<-p3 + geom_text(aes(label = mylabels))
ggpubr::ggarrange(p1)
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

data<-read.csv("Palicourea_median_climate_data_with_elev_pollination_final.csv")
#mean(data$elev)

pollination<-data$Principal_pollinator
species<-data$species
pollination<- setNames(pollination, species)

elevation<-data$elev
elevation<- setNames(elevation, species)
astral01 <- read.tree("../Resubmission/rev_dendrogram.tre")

tips2delete<-c("Pal_guianensis2","Pal_demissa2","Pal_angustifolia2","Pal_obliquinervia",
               "Pal_wolffiae","Pal_cauligera","Pal_vesiculifera","Pal_alagoana",
               "Pal_laxivenulosa","Pal_angustifolia1","Pal_wolffiae","Pal_cauligera")
tree <- ape::drop.tip(astral01, tips2delete)
plot(tree)
new_names <- read.csv("~/Repos/Palicourea/Climatic_niche_overlap/pal_name_updates.csv", header = T)
match_indices <- match(tree$tip.label, new_names$old_label)
# Replace names in tree$tip.label vector
tree$tip.label <- ifelse(is.na(match_indices), tree$tip.label, new_names$new_label[match_indices])
plot(tree)

#Check that tree and data have same species
tree_species <- tree$tip.label
data_species <- unique(data$species)               


data_not_in_tree <- setdiff(data_species, tree_species)
length(data_not_in_tree)
head(data_not_in_tree)

tree_not_in_data <- setdiff(tree_species, data_species)
length(tree_not_in_data)
head(tree_not_in_data)
#convert Species column into row names

phylANOVA(tree, pollination, elevation)

p <- ggplot(data, aes(x=Principal_pollinator, y=elev, color=Principal_pollinator)) + # fill=name allow to automatically dedicate a color for each group
  geom_boxplot()
p

library("car")

levene_test = leveneTest(elev ~ Principal_pollinator, data)

print(levene_test)
#phylANOVA(tree, pollination, bio12)

##reading data for proportion of pollinators at given elevation
df<-read.csv("Palicourea_raw_climate_data_with_elev_pollinators_prop.csv")


bin_width <- 50

result <- df %>%
  filter(!is.na(elev), !is.na(species), !is.na(pollinator)) %>%
  mutate(
    # Bin index
    bin_id = case_when(
      elev <= 50 ~ 0,
      TRUE ~ floor((elev - 1) / bin_width)
    ),
    
    # Bin boundaries
    bin_start = case_when(
      bin_id == 0 ~ 0,
      TRUE ~ bin_id * bin_width + 1
    ),
    bin_end = (bin_id + 1) * bin_width,
    
    elev_bin = paste0(bin_start, "-", bin_end),
    
    is_hummingbird = grepl("hummingbird", pollinator, ignore.case = TRUE)
  ) %>%
  distinct(species, bin_id, elev_bin, bin_start, bin_end, is_hummingbird) %>%
  group_by(bin_id, elev_bin, bin_start, bin_end) %>%
  summarise(
    n_species = n_distinct(species),
    n_hummingbird_species = n_distinct(species[is_hummingbird]),
    prop_hummingbird = n_hummingbird_species / n_species,
    .groups = "drop"
  ) %>%
  arrange(bin_id)

write_csv(result, "hummingbird_proportion_by_elevation_bin.csv")

proportions<-read_csv("hummingbird_proportion_by_elevation_bin.csv")
reggression<-lm(prop_hummingbird ~ bin_end, data = proportions)
summary(reggression)

ggplot(result, aes(x = bin_end, y = prop_hummingbird)) +
  geom_point(size = 2, color = "grey70") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(
    x = "Elevation (m)",
    y = "Proportion (hummingbird pollinated / all species)"
  ) +
  theme_minimal()
############################
############################
###########################






###########################################################################
##Testing for the effect of spp. with broad elevation ranges on phylANOVA##
###########################################################################

niche_data_complete<-read.csv("Palicourea_raw_climate_data_with_elev_pollinators_prop.csv")

#Estimate standard deviation (SD) of elevation for each species
sd_niche_values<-aggregate(niche_data_complete[2:3], niche_data_complete[1], sd)
write.csv(sd_niche_values, "Palicourea_sd_climate_data_with_elev_pollinators_prop.csv", row.names=F)

#Remove species with SD above the average SD across species
mean(sd_niche_values$elev) #438.8847

broad_elev_removed = sd_niche_values[sd_niche_values$elev > 438.8847,]

new_names <- read.csv("pal_name_updates.csv", header = T)
match_indices <- match(tree$tip.label, new_names$old_label)
# Replace names in tree$tip.label vector
tree$tip.label <- ifelse(is.na(match_indices), tree$tip.label, new_names$new_label[match_indices])

spp2remove<-broad_elev_removed$species

new_data <- subset(data, !(species %in% spp2remove))

new_pollination<-new_data$Principal_pollinator
new_species<-new_data$species
new_pollination<- setNames(new_pollination, new_species)

new_elevation<-new_data$elev
new_elevation<- setNames(new_elevation, new_species)

new_tree <- ape::drop.tip(tree, spp2remove)

phylANOVA(new_tree, new_pollination, new_elevation)

new_p <- ggplot(new_data, aes(x=Principal_pollinator, y=elev, color=Principal_pollinator)) + # fill=name allow to automatically dedicate a color for each group
  geom_boxplot()
new_p

library("car")

new_levene_test = leveneTest(elev ~ Principal_pollinator, new_data)

print(new_levene_test)


