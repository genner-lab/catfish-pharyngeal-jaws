library(tidyverse)
library(ggtree)
library(ape)
library(tidytree)
library(phytools)
library(ggtreeExtra)
library(treeio)
library(here)
library(readr)
library(dplyr)
setwd("C:/Users/sl22263/OneDrive - University of Bristol/Documents/Data to be landmarked/phylogeny")

# read ml tree, bootstrap trees, species list
tr <- ape::read.tree(("bootstrapped.raxml.support"))
bs <- ape::read.tree(("concatenated.aligned.trimmed.fasta.raxml.bootstraps"))
keep <- readr::read_csv(("Species.csv"),show_col_types=FALSE) |> dplyr::pull(Name)

# tree pruning and rooting function
tree_prune <- function(tr,tips,root) {
  tr.d <- ape::keep.tip(tr,tips)
  og <- grep(root,tr.d$tip.label,value=TRUE)
  tr.r <- ape::root(tr.d,outgroup=og,resolve.root=TRUE)
  return(tr.r)
}

# prune all bootstrap trees
bs.sub.root <- lapply(bs,function(x) tree_prune(tr=x,tips=keep,root="Astroblepus"))

# prune ml tree
tr.sub.root <- tree_prune(tr=tr,tips=keep,root="Astroblepus")

# recalculate bootstraps
boot <- ape::prop.clades(tr.sub.root,bs.sub.root,rooted=TRUE)/10

# transfer bootstraps to tree and remove root support
tr.sub.root$node.label <- boot
tr.sub.root$node.label[1] <- ""

# write out
ape::write.tree(tr.sub.root,here::here("bootstrapped.raxml.support.recalc.nwk"))

tree <- read.tree("bootstrapped.raxml.support.recalc.nwk")
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

traits <- namecorrected %<+% jaw 
  
####ancestral states
jaw <- read.csv("pharyngealjawtypenodiet.csv",row.names=1)
jaw<-as.matrix(jaw)[,1]
jaw
namecorrected$node.label <- NULL
library(geiger)
fitDiscrete(namecorrected, jaw,
            model = "SYM")
#####aic summary:
#log-likelihood = -35.511902
#AIC = 101.023804
#AICc = 111.023804
#free parameters = 15
fitDiscrete(namecorrected, jaw,
            model = "ER")
###AIC summary
#log-likelihood = -43.963604
#AIC = 89.927209
#AICc = 89.991725
#free parameters = 1

fitER<-ace(jaw,namecorrected,model="ER",type="discrete", CI = TRUE)
fitER
  
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
final <- ggtree(rotated) + geom_tiplab() +
  geom_point2(aes(colour = "black", label = label, size = 0.8, subset = !is.na(as.numeric(label)) & as.numeric(label) >70)) +
  xlim(NA, 0.8) 
plot(final)
########################################
final <- ggtree(namecorrected) + xlim(NA, 0.8) 
plot(final)
#flip to Planiloriaria at top. 
rotated <- flip(final, 85, 70)
plot(rotated)
final <- ggtree(rotated) + geom_tiplab() +
  geom_point2(aes(colour = "black", label = label, size = 0.8, subset = !is.na(as.numeric(label)) & as.numeric(label) >70)) +
  xlim(NA, 0.8) 
plot(final)

#mapping traits
jaw <- read.csv("pharyngealjawtypenodiet.csv")
jaw$Jaw <- as.factor(jaw$Jaw)
final <- rotated %<+% jaw + geom_tippoint(aes(colour = Jaw, label = NULL), position = position_nudge(0.007), size = 2.0) +
  scale_colour_manual(values = c('Labroid-like' = "cyan3", 'Elongated, elongate teeth'="purple", 'Bladed, no teeth'="green", 'Bladed, small teeth'="red", 'Comb-shaped'="orange", 'Carp-like'="grey")) +
  guides(colour = "none") +
  theme(legend.position = "top left") 

plot(final)
x<-getStates(rotated1,"jaw")

#basic tree 
ggtree(namecorrected) +
  geom_nodelab(aes(label=label, check_overlap = TRUE), nudge_x=0.008, size = 3.0) +
  geom_tiplab(geom = "text", aes(fontface = "italic")) +
  xlim(NA, 0.8) 




