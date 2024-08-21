library(tidyverse)
library(ggtree)
library(ape)
library(tidytree)
library(phytools)
library(ggtreeExtra)
library(treeio)
setwd("C:/Users/sl22263/OneDrive - University of Bristol/Documents/Data to be landmarked/phylogeny")
tree <- read.tree("bootstrapped.raxml.support")
ggplot(tree) + geom_tree() + theme_tree() + geom_nodelab(aes(label=node))
ggtree(tree)
# create the basic plot
p <- ggtree(tree)
# add node points
p + geom_nodepoint()
# add tip points
p + geom_tippoint()
# Label the tips
p + geom_tiplab()
#rooting the tree
outgroup <- MRCA(as.phylo(p), "Astroblepus_caquetae") 
rooted_tree <- root(as.phylo(tree), outgroup)
#prune nodes
species<-read.csv("Species.csv", header = T)
species <- c("Acanthicus_hystrix", "Ancistrus_brevipinnis", "Ancistrus_ranunculus", "Aphanotorulus_emarginatus", "Astroblepus_caquetae", "Baryancistrus_chrysolomus", "Baryancistrus_demantoides", "Chaetostoma_dermorhynchus", "Cordylancistrus_platycephalus",
             "Crossoloricaria_variegata", "Cteniloricaria_platystoma", "Dekeyseria_scaphirhynchus", "Farlowella_acus", "Harttia_loricariformis", "Hemiancistrus_medians", "Hemiodontichthys_acipenserinus", "Hisonotus_notatus",
             "Hypancistrus_debilittera", "Hypancistrus_furunculus", "Micracanthicus_vandragti", "Hypoptopoma_guianense", "Hypoptopoma_incognitum", "Hypostomus_cochliodon", "Hypostomus_commersoni", "Hypostomus_laplatae", "Hypostomus_plecostomus",
             "Hypostomus_regani", "Hypostomus_watwata", "Isorineloricaria_spinosissima", "Lamontichthys_filamentosus", "Lasiancistrus_caucanus", "Leporacanthicus_heterodon", "Limatulichthys_griseus", "Lithogenes_villosus",
             "Lithoxus_lithoides", "Loricaria_cataphracta", "Loricariichthys_maculatus", "Neoplecostomus_franciscoensis", "Otocinclus_vittatus", "Otothyris_juquiae", "Oxyropsis_carinata",
             "Panaqolus_maccus", "Panaque_armbrusteri", "Parancistrus_nudiventris", "Pareiorhaphis_calmoni", "Parotocinclus_jumbo", "Peckoltia_braueri", "Peckoltia_sabaji",
             "Planiloricaria_cryptodon", "Pseudacanthicus_sp.Xingu_L025", "Pseudancistrus_barbatus", "Pseudohemiodon_apithanos", "Pseudoloricaria_laeviuscula", "Pseudorinelepis_genibarbis", "Pterygoplichthys_gibbiceps", "Parotocinclus_eppleyi",
             "Rhinelepis_aspera", "Rineloricaria_lanceolata", "Scobinancistrus_aureatus", "Spatuloricaria_cf.evansii_LBP", "Spectracanthicus_immaculatus", "Spectracanthicus_zuanoni", "Sturisoma_robustum", "Sturisomatichthys_tamanae")
pruned.tree <- keep.tip(rooted_tree, species)
ggtree(pruned.tree)
#changed node names
ID <- read.csv("Species.csv", header = T)
ID1 <-as.data.frame(ID)
namecorrected <-rename_taxa(pruned.tree, ID1, key = 1, value = 2)
plot(namecorrected)
#final tree
final <- ggtree(namecorrected) + geom_tiplab(geom = "text", aes(fontface = "italic")) +
  geom_nodelab(aes(label = label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70
  )) 
plot(namecorrected)   
#no bootstraps
final <- ggtree(namecorrected) + geom_tiplab(geom = "text", aes(fontface = "italic"), nudge_x=0.01) +
  geom_point2(aes(colour = "black", label = label, size = 0.8, subset = !is.na(as.numeric(label)) & as.numeric(label) >70)) +
  xlim(NA, 0.8) 
plot(final)
#flip to Planiloriaria at top. 
rotated <- flip(final, 85, 70)
plot(rotated)

#mapping traits
jaw <- read.csv("pharyngealjawtypenodiet.csv")
jaw$Jaw <- as.factor(jaw$Jaw)
final <- rotated %<+% jaw + geom_tippoint(aes(colour = Jaw, label = NULL), position = position_nudge(0.007), size = 2.0) +
  scale_colour_manual(values = c('Labroid-like' = "cyan3", 'Elongated, elongate teeth'="purple", 'Bladed, no teeth'="green", 'Bladed, small teeth'="red", 'Comb-shaped'="orange", 'Carp-like'="grey")) +
  guides(colour = "none") +
  theme(legend.position = "top left") 

plot(final)

#basic tree 
ggtree(namecorrected) + 
  geom_nodelab(aes(label=label)) +
  geom_tiplab(geom = "text", aes(fontface = "italic"), nudge_x=0.01)




