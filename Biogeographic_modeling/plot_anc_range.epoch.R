source("/Applications/RevBayes_OSX_v1.0.13/scripts/plot_anc_range.util.R")

setwd("~/Desktop/Palicourea/RevBayes/simple/")
# file names
fp = "~/Desktop/Palicourea/RevBayes/simple/" # edit to provide an absolute filepath
plot_fn = paste(fp, "simple.range.pdf",sep="")
mcc_fn = paste(fp, "rev_DTE.MCC.tre", sep="")
tree_fn = paste(fp, "rev_DTE.ase.tre", sep="")
label_fn = paste(fp, "rev_DTE.state_labels.txt", sep="")
color_fn = paste(fp, "pali_range_colors.txt", sep="")


# get state labels and state colors
states = make_states(label_fn, color_fn, fp=fp)
state_labels = states$state_labels


# process the ancestral states
ase <- processAncStates(tree_fn,
                        # Specify state labels.
                        # These numbers correspond to
                        # your input data file.
                        state_labels = state_labels)


ncol <- length(state_labels)
# C (Central America), L (Clow inter Andean), A (High Andes), E (Eastern Amazon), O (Orinoco), F(Atlantic Forest)
colors_main <- c("#4E598C", "#FB9B2D", "#4ECDC4", "#568259","#FF6B6B","#97CC04")

colors_combined <- colorRampPalette(c("#FDA8EF","cyan", "red","royalblue",
                                      "aquamarine", "#73C1FC", "darkviolet"))(128-7)


colors <- c(colors_main, colors_combined)
names(colors) <- state_labels

pp_map=plotAncStatesMAP(t = ase, node_color = colors, tip_labels_size = 2.7,
                        tip_labels_italics = T, tip_labels_states = T,
                        node_labels_as = NULL, node_labels_centered = T,
                        node_labels_offset = 0, tip_labels_states_offset = 1,
                        tip_labels_states_size = 2, tip_states=F,
                        cladogenetic = TRUE, tip_labels_offset = 2, timeline = T, geo_units="epochs", tip_age_bars=TRUE, node_age_bars=TRUE, age_bars_colored_by="posterior", label_sampled_ancs=TRUE)+
  ggplot2::theme(legend.position = c(0.08, 0.67),
                 legend.key = element_blank(), legend.title = element_blank(),
                 legend.text = element_text(size = 11),
                 legend.background = element_rect(fill="transparent"))+
  
  # andean uplift peaks 
  ggplot2::geom_vline(xintercept = -23, linetype = "dotted") +
  ggplot2::annotate(geom = "text", 
                    x = -23, y = 97, 
                    label = "A1",
                    hjust = 1, 
                    size = 8)+
  ggplot2::geom_vline(xintercept = -12, linetype = "dotted") +
  ggplot2::annotate(geom = "text", 
                    x = -12, y = 97, 
                    label = "A2",
                    hjust = 1, 
                    size = 8)+
  ggplot2::geom_vline(xintercept = -4.5, linetype = "dotted") +
  ggplot2::annotate(geom = "text", 
                    x = -4.5, y = 97, 
                    label = "A3",
                    hjust = 1, 
                    size = 8)
#save
ggsave(file=plot_fn, plot=pp_map, device="pdf", height=15, width=12, useDingbats=F)

##Plotting divergence times with 95% HPD intervals. Modofied from David Černý (https://davidcerny.github.io/post/plotting_beast/)

library(phytools)
library(strap)
library(phyloch)

# Get shaded bars for the HPD intervals. Credit:
# http://blog.phytools.org/2017/03/error-bars-on-divergence-times-on.html

annot_tree <- phyloch::read.beast("~/Desktop/Palicourea/Manuscript/calibration/standard_Palicourea_ORC_versioncorrected_MCC.tre")
tips2delete<-c("Car_guianensis_Gonzalez2158","Carapichea_ipecacuanha_Croat15117",
               "Eum_boliviana_Campbell22035","Not_epiphytica_Neill15737",
               "Not_uliginosa_Stevens37138","Psy_carthagenensis_Araujo2124",
               "Psy_grandis_Taylor11745","Psy_guianensis_Merello1711",
               "Psy_horizontalis_Stevens32733","Psy_jinotegensis_Stevens33549",
               "Psy_limonensis_Stevens31580","Psy_marginata_Stevens32781",
               "Psy_nervosa_Stevens32362","Psy_panamensis_Stevens32285","Psy_subsessilis_Stevens31494","Rud_cornifolia_deGracias818")

if (is.null(annot_tree$`CAheight_95%_HPD_MIN`)) {
  annot_tree$min_ages <- annot_tree$`height_95%_HPD_MIN`
  annot_tree$max_ages <- annot_tree$`height_95%_HPD_MAX`
} else {
  annot_tree$min_ages <- annot_tree$`CAheight_95%_HPD_MIN`
  annot_tree$max_ages <- annot_tree$`CAheight_95%_HPD_MAX`
}

annot_tree$root.time <- max(nodeHeights(annot_tree)) + 0.0
time_tree=geoscalePhylo(ladderize(annot_tree, right = F), x.lim = c(0, 45), cex.tip = 0.7, cex.age = 1.3, cex.ts = 1)

T1 <- get("last_plot.phylo", envir = .PlotPhyloEnv)

for(i in (Ntip(annot_tree) + 1):(annot_tree$Nnode + Ntip(annot_tree))) {
  lines(x = c(T1$root.time - annot_tree$min_ages[i - Ntip(annot_tree)],
              T1$root.time - annot_tree$max_ages[i - Ntip(annot_tree)]),
        y = rep(T1$yy[i], 2), lwd = 4, lend = 0,
        col = make.transparent("blue", 0.3))
}

#misc
# plot the ancestral states
#pp=plotAncStatesPie(t=ase_Pali,
#                      include_start_states=T,
#                      summary_statistic="PieRange",
#                      state_labels=state_labels,
#                      state_colors=colors,
#                      tip_label_size=2.5,
#                      tip_label_offset=0.1,
#                      node_label_size=0,
#                      shoulder_label_size=0,
#                      show_posterior_legend=T,
#                      tip_pie_diameter=0.5,
#                      node_pie_diameter=2.0,
#                      pie_nudge_x=0.03,
#                      pie_nudge_y=0.16,
#                      alpha=1,
#                      timeline=T)
# Move the legend
#  theme(legend.position = c(0.1, 0.75))  
# get plot dimensions
#x_phy = max(pp$data$x)       # get height of tree
#x_label = 10                # choose space for tip labels
#x_start = 40                  # choose starting age (greater than x_phy)
#x0 = -(x_start - x_phy)      # determine starting pos for xlim
#x1 = x_phy + x_label         # determine ending pos for xlim

#pdf(paste0("plot_pies.pdf"),height = 20, width = 20)
#pp=plotAncStatesPie(t = ase, pie_colors = colors, 
#                      tip_pie_nudge_x = 1, tip_labels_states_offset = 0.5,
#                      tip_labels_size = 3, node_pie_size = 1,
#                      tip_pie_size = 0.4, tip_labels_states = T,
#                      cladogenetic = TRUE, tip_labels_offset = 2, 
#                      timeline = T) +
