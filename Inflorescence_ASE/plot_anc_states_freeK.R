################################################################################
#
# Plot ancestral state estimates inferred using an independent rates models of character evolution.
#
#
# authors: Sebastian Hoehna
#
################################################################################

library(RevGadgets)
library(ggplot2)

CHARACTER <- "palicourea.inflorescence_type"
NUM_STATES <- 4

STATE_LABELS <- c("0" = "psy_open", "1" = "psy_congested", "2" = "Pali", "3" = "Cephaelis")

tree_file <- paste0(CHARACTER,"_ase_freeK.tree")

# process the ancestral states
ase <- processAncStates(tree_file,
                        # Specify state labels.
                        # These numbers correspond to
                        # your input data file.
                        state_labels = STATE_LABELS)

# produce the plot object, showing MAP states at nodes.
# color corresponds to state, size to the state's posterior probability
p <- plotAncStatesMAP(t = ase,
                      tree_layout = "rect",
                      tip_labels_size = 1) +
     # modify legend location using ggplot2
     theme(legend.position = c(0.92,0.81))

ggsave(paste0(CHARACTER,"_ASE_FreeK_MAP.pdf"), p, width = 11, height = 9)




p <- plotAncStatesPie(t = ase,
                      tree_layout = "rect",
                      tip_labels_size = 1) +
     # modify legend location using ggplot2
     theme(legend.position = c(0.92,0.81))

ggsave(paste0(CHARACTER,"_ASE_FreeK_Pie.pdf"), p, width = 11, height = 9)
