library(ape)
tree <- read.tree("ASTRAL_810_allexons-BS10.nwk")
tips_to_drop <- c("COU_hexandra_Callejas4587","Pal_elata_Stevens36104","Pal_correae_Dwyer1968",
                  "Pal_tetragona_Davidse36909","Pal_suerrensis_McPherson15864",
                  "Pal_acuminata_Suazo4616","Pal_triphylla_Nee46793","Pal_triphylla_Killeen6548",
                  "Pal_rigida_Subieta322", "Pal_deflexa_Meave1172","Pal_pubescens_Stevens21197",
                  "Pal_brachiata_Taylor11718", "HOF_phoenicopoda_Wendt3392","HAM_patens_LB54",
                  "COS_grandiflora_Tyson904","HIL_parasitica_Jiminez2189","BAL_stormiae_Vasquez1081")

tree <- ape::drop.tip(tree, tips_to_drop)
#tree <- multi2di(tree) # remove multifurcations
#tree$edge.length <- pmax(tree$edge.length,1/365)
#di2multi(phy, tol = 1e-08)
rooted_tree<-root(tree, outgroup=c("Psy_jinotegensis_Stevens33549","Psy_nervosa_Stevens32362",
                                   "Psy_grandis_Taylor11745","Psy_panamensis_Stevens32285",
                                   "Psy_limonensis_Stevens31580","Psy_subsessilis_Stevens31494",
                                   "Psy_horizontalis_Stevens32733","Psy_marginata_Stevens32781",
                                   "Psy_guianensis_Merello1711","Psy_carthagenensis_Araujo2124"), resolve.root=T)

plot(rooted_tree)


#tips_to_drop <- c("BAL_stormae_Vasquez1081","COS_grandiflora_Tyson904","HAM_patens_LB54","HOF_phoenicopoda_Wendt3392","COU_hexandra_Callejas4587")

#write.tree(tree,file = "ASTRAL_188singlecopy-BS10_rooted_outgroup.tre")

# date tree with PL 

mycalibration <- data.frame(
  node = c( # Luzuriaga stem age fossil calibration
    #getMRCA(tree, tip = c("HOF_phoenicopoda_Wendt3392","HAM_patens_LB54","COS_grandiflora_Tyson904","HIL_parasitica_Jiminez2189","BAL_stormae_Vasquez1081")),
    # Psychotrieae crown age secondary calibration 
    getMRCA(tree, tip=c("Psy_horizontalis_Stevens32733", "Pal_domingensis_Axelrod1028"))),
    #getMRCA(tree, tip = c("Bomarea_salsilla_Chile_Ackerman545",
    #                      "Bomarea_multiflora_CultivatedinCAfromCol_Greenhouse")),
    # Alstroemeria crown age secondary calibration 
    #getMRCA(tree, tip = c("Alstroemeria_apertiflora_Hatschbach17552", 
    #                      "Alstroemeria_revoluta_Watson6608")),
    # root age secondary calibration 
    #getMRCA(tree, tip = c("Luzuriaga_polyphylla_CultivatedinCA_UCBG90_2401", 
    #                      "Bomarea_salsilla_Chile_Ackerman545"))),
  age.min = 20,
  age.max = 46.9,
  soft.bounds = TRUE)
tree_dated <- chronos(rooted_tree, lambda = 2, model = "relaxed",
                        calibration = mycalibration,
                        control = chronos.control() )
lad<-ladderize(tree_dated,right=F)
plot(lad)
write.tree(tree_dated, file = "chronos_ASTRAL_all_exons_810-BS10.tre")