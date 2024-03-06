#Script modified by Laymon Ball. Original script created by Shuiyin Liu at 15:43, FRI, 2023-11-17 (https://github.com/ShuiyinLIU/QS_visualization/blob/master/plot_QC_ggtree.R)

setwd("~/Desktop/Pal_QSC//")
list.files()

library(ggtree)
library(treeio)
library(ggplot2)
library(ape)

qc <- read.tree("RESULT.labeled.tre.qc")
qd <- read.tree("RESULT.labeled.tre.qd")
qi <- read.tree("RESULT.labeled.tre.qi")


# process node labels of above three labeled trees
# qc tree
tree <- qc
tree$node.label <- gsub("qc=","",tree$node.label)
#View(tree$node.label)
write.tree(tree,"tree_qc.tre")

# qd tree
tree <- qd
tree$node.label <- gsub("qd=","",tree$node.label)
#View(tree$node.label)
write.tree(tree,"tree_qd.tre")

# qi tree
tree <- qi
tree$node.label <- gsub("qi=","",tree$node.label)
#View(tree$node.label)
write.tree(tree,"tree_qi.tre")


# read 3 modified tree files for QC/QD/QI
tree_qc <- read.newick("tree_qc.tre", node.label='support')
tree_qd <- read.newick("tree_qd.tre", node.label='support')
tree_qi <- read.newick("tree_qi.tre", node.label='support')


# add a customized label for internode or inter-branch, i.e., qc/qd/qI
score_raw = paste(tree_qc@data$support,"/",tree_qd@data$support,"/",tree_qi@data$support,sep="")
score = gsub("NA/NA/NA","",score_raw)
score = gsub("NA","-",score)
#View(score)


# set labeled QC tree as the main plot tree
tree <- tree_qc
tree@data$score <- score


#####################################################
# Partitioning quartet concordance. QC values were divided into four categories and this
# information was used to color circle points. 

# drop the internodes without QC vaule
root <- tree@data$node[is.na(tree@data$support)]

#pdf(file="00.treeQC_rect_circ.pdf", width = 16, height = 18)
pdf(file="pies.pdf", width = 16, height = 18)

#view tree with node numbers
ggtree(tree) + geom_text(aes(label = node))



####(1)color circle points####
ggtree(tree, size=0.5, branch.length = "none") +
  geom_tiplab(size=3) + xlim(0, 22) +
  geom_nodepoint(aes(subset=!isTip & node != root, fill=cut(support, c(1, 0.2, 0, -0.05, -1))),
                 shape=21, size=4) +
  theme_tree(legend.position=c(0.15, 0.88)) +
  scale_fill_manual(values=c("#117733", "#C0DA71", "#FFCC99","#EA5E00"),
                    guide="legend", name="Quartet Concordance(QC)",
                    breaks=c("(0.2,1]","(0,0.2]","(-0.05,0]","(-1,-0.05]"),
                    labels=expression(QC>0.2, 0 < QC * " <= 0.2", -0.05 < QC * " <= 0", QC <= -0.05))

####(2)color branches and add pie charts####
tree2 <- tree #read.newick("tree_qc.tre")
tree2@phylo$data <- as.data.frame(read.csv("QS_results_combined.csv"))

qd <- as.numeric(unlist(read.csv("QS_results_combined.csv")[9]))
qd <- round(qd , digits = 2)
tree2@data$support_qd <- qd

#p <- ggtree(tree2)#, branch.length = "none", size = 0.5) + geom_text2(aes(label = node)) + geom_tiplab() + xlim(0,35)
pies <- nodepie(tree2@phylo$data, cols = 5:7, alpha = 0.9, color = c(count0_prop='#332288', count1_prop='#56B4E9', count2_prop='#D81B60'))

tree3 <- ggtree(tree2, aes(color=cut(support, c(1, 0.2, 0, -0.05, -1))), size=0.75, branch.length = "none") +
  geom_tiplab(size=2) + xlim(0, 22) +
  #geom_text2(aes(label = support_qd, x = branch, fontface = 2), size = 3, vjust = -.5, color = "black") +
  theme_tree(legend.position=c(0.15, 0.88)) +
  guides(color = guide_legend(override.aes = aes(label = ""))) +
  scale_colour_manual(values=c("#117733", "#C0DA71", "#FFCC99","#EA5E00"),
                      breaks=c("(0.2,1]","(0,0.2]","(-0.05,0]","(-1,-0.05]"),
                      na.translate=T, na.value="black",
                      guide="legend", name="Quartet Concordance(QC)",
                      labels=expression(QC>0.2, 0 < QC * " <= 0.2", -0.05 < QC * " <= 0", QC <= -0.05)) 
#tree3

#black branches
tree4 <- ggtree(tree2, size=1, branch.length = "none") +
  geom_tiplab(size=3) + xlim(0, 22) +
  #geom_text2(aes(label = support_qd, x = branch, fontface = 2), size = 3, vjust = -.5, color = "black") +
  theme_tree(legend.position=c(0.15, 0.88)) 

inset(tree3, pies, width = 0.035, height = 0.035, hjust = 0.01, vjust = 0.15) #adjust width and height to change size of pie charts 


####(3)color circle points and label each internode with QC/QD/QI####
ggtree(tree, size=0.5, branch.length = "none") +
  geom_tiplab(size=3) + xlim(0, 22) +
  geom_nodepoint(aes(subset=!isTip & node != root, fill=cut(support, c(1, 0.2, 0, -0.05, -1))),
                 shape=21, size=4) +
  theme_tree(legend.position=c(0.15, 0.88)) +
  scale_fill_manual(values=c("#117733", "#C0DA71", "#FFCC99","#EA5E00"),
                    guide="legend", name="Quartet Concordance(QC)",
                    breaks=c("(0.2,1]","(0,0.2]","(-0.05,0]","(-1,-0.05]"),
                    labels=expression(QC>0.2, 0 < QC * " <= 0.2", -0.05 < QC * " <= 0", QC <= -0.05))+
  geom_text(aes(label=score, x=branch), size=3, vjust=-.5)

dev.off()





