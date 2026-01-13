library(geomorph) #
library(ape) #
library(ggplot2) #
library(MASS) # 
library(ggpubr) #
library(dplyr)
library(MASS)
library(caret)
library(tidyr)
library(Morpho)
setwd("C:/Users/sl22263/OneDrive - University of Bristol/Documents/Data to be landmarked/CT scan images/Loricariidae Samples/analyses")
tree <- read.nexus(file = "phylogenyLoricariidae2.nex")
speciephy <- c("Ancistrus brevipinnis", "Ancistrus ranunculus", "Aphanotorulus emarginatus", "Baryancistrus chrysolomus", "Baryancistrus demantoides", "Chaetostoma dermorhynchus",
               "Dekeyseria scaphirhynchus", "Farlowella acus", "Harttia loricariformis", "Hisonotus notatus",
               "Hypancistrus debilittera", "Hypancistrus inspector", "Hypoptopoma psilogaster", "Hypostomus cochliodon", "Hypostomus commersoni", "Hypostomus laplatae", "Hypostomus plecostomus",
               "Hypostomus regani", "Hypostomus watwata", "Lamontichthys filamentosus", "Lasiancistrus mayoloi", "Leporacanthicus heterodon", "Limatulichthys griseus",
               "Lithoxus lithoides", "Loricaria cataphracta", "Loricariichthys maculatus", "Neoplecostomus granosus", "Otocinclus vittatus",
               "Panaqolus maccus", "Panaque armbrusteri", "Parancistrus aurantiacus", "Pareiorhaphis cameroni", "Peckoltia braueri",
               "Planiloricaria cryptodon", "Pseudacanthicus pirarara", "Pseudohemiodon apithanos", "Pseudoloricaria laeviuscula", "Pterygoplichthys gibbiceps", "Rhinotocinclus eppleyi",
               "Rhinelepis strigosa", "Rineloricaria lanceolata", "Scobinancistrus aureatus", "Spatuloricaria evansii", "Spectracanthicus zuanoni", "Sturisoma robustum")
tree <- keep.tip(tree, speciephy)

###identifiers
ID <- read.csv("IDpharyphy no outgroup.csv", header = T)
pharycat <-factor(ID$oral.jaw.1)
ID <- ID[!(ID$Niche %in% "Unknown"),]
##pharyngeal jaw shape
phary <- factor(ID$pharyngealjaw)
##oral jaw shape
oralinfo <- factor(ID$oral.jaw.1)

###Diet identifiers 
ID <- read.csv("IDdietfinalcat.csv", header = T)
categories <-factor(ID$cat)
Species <-factor(ID$name)

ID1 <- subset(ID, select = -c(Species, cat, special, Location) )
library(vegan)
sol <- metaMDS(ID1, distance = "euclidean")
sol
plot(sol, type="t")
MDS <- as.data.frame(sol$points)
ID <- cbind(ID, MDS)

plotmds <- ggscatter(ID, x = "MDS1", y = "MDS2", color = "cat", label = "Species", font.label = c("italic", "black"),
                     ellipse = TRUE, ellipse.type = "convex", ellipse.alpha = 0.08, mean.point = FALSE, point = "FALSE", star.plot = FALSE, legend.title = "Dietary Niche") + 
  theme_classic() +
  geom_point() +
  labs(x ="PC1 (X% of variation)", y = "PC2 (X% of variation)") +
  theme(legend.position = "right") 
plotmds

#####MDS corrected identifiers
ID <- read.csv("IDdietfmds.csv", header = T)
categories <-factor(ID$cat)
Species <-factor(ID$name)

ID$phary <- factor(phary)
##oral jaw shape
ID$oralinfo <- factor(oralinfo)

ID1 <- subset(ID, select = -c(Species, cat, special, Location, phary, oralinfo) )
sol <- metaMDS(ID1, distance = "euclidean")
sol
plot(sol, type="t")
MDS <- as.data.frame(sol$points)
ID <- cbind(ID, MDS)

plotmds <- ggscatter(ID, x = "MDS1", y = "MDS2", color = "cat", font.label = c("italic", "black"),
                     ellipse = TRUE, ellipse.type = "convex", ellipse.alpha = 0.08, mean.point = FALSE, point = "FALSE", star.plot = FALSE, legend.title = "Dietary Niche",
                     palette = c('Detritivore' = "blue", 'Algivore' = "limegreen", 'Insectivore' = "orange")) + 
  theme_classic() +
  geom_point() +
  labs(x ="NMDS1", y = "NMDS2") +
  theme(legend.position = "right") 
plotmds

IDalgae <- read.csv("IDalgae.csv", header = T)
IDDetritus <- read.csv("IDDetritus.csv", header = T)
IDcarnivore <- read.csv("IDcarnivore.csv", header = T)


d <- ggbarplot(IDDetritus, "Category", "Mean", color = "Black", fill = "Black", 
               legend = "Null") +
  scale_y_continuous(limits=c(0, 21), breaks = seq(0, 30, by = 1),expand = expansion(mult = c(0, 0.05))) +
  theme(axis.text.x = element_blank()) +
  labs(x =NULL, y = NULL) +
  coord_flip()
d

a <- ggbarplot(IDalgae, "Category", "Mean", color = "Black", fill = "Black", 
               legend = "Null") +
  theme(axis.text.x = element_blank()) +
  labs(x =NULL, y = NULL) +
  scale_y_continuous(limits=c(0, 21), breaks = seq(0, 1000, by = 1), expand = expansion(mult = c(0, 0.05))) +
  coord_flip()
a
c <- ggbarplot(IDcarnivore, "Category", "Mean", color = "Black", fill = "Black", 
               legend = "Null") +
  theme(axis.text.x = element_blank()) +
  labs(x =NULL, y = NULL) +
  scale_y_continuous(limits=c(0, 21), breaks = seq(0, 1000, by = 1), expand = expansion(mult = c(0, 0.05))) +
  coord_flip()
c
dietcomponents <- ggarrange(d, a, c, nrow = 3, align = "hv")
dietcomponents



#####oral jaws#########

#######comparing oral jaw classifiers to oral jaw shape 
treeout <- read.nexus(file = "phylogenyLoricariidae2.nex")
library(tidytree)
speciephynooutgroup <- c("Acanthicus hystrix", "Ancistrus brevipinnis", "Ancistrus ranunculus", "Aphanotorulus emarginatus", "Baryancistrus chrysolomus", "Baryancistrus demantoides", "Chaetostoma dermorhynchus", "Cordylancistrus setosus",
                         "Crossoloricaria variegata", "Cteniloricaria platystoma", "Dekeyseria scaphirhynchus", "Farlowella acus", "Harttia loricariformis", "Hemiancistrus medians", "Hemiodontichthys acipenserinus", "Hisonotus notatus",
                         "Hypancistrus debilittera", "Hypancistrus inspector", "Hypancistrus vandragti", "Niobichthys ferrarisi", "Hypoptopoma psilogaster", "Hypostomus cochliodon", "Hypostomus commersoni", "Hypostomus laplatae", "Hypostomus plecostomus",
                         "Hypostomus regani", "Hypostomus watwata", "Isorineloricaria spinosissima", "Lamontichthys filamentosus", "Lasiancistrus mayoloi", "Leporacanthicus heterodon", "Limatulichthys griseus", "Lithogenes villosus",
                         "Lithoxus lithoides", "Loricaria cataphracta", "Loricariichthys maculatus", "Neoplecostomus granosus", "Otocinclus vittatus", "Otothyris sp.", "Oxyropsis carinata",
                         "Panaqolus maccus", "Panaque armbrusteri", "Parancistrus aurantiacus", "Pareiorhaphis cameroni", "Parotocinclus jumbo", "Peckoltia braueri", "Peckoltia sabaji",
                         "Planiloricaria cryptodon", "Pseudacanthicus pirarara", "Pseudancistrus barbatus", "Pseudohemiodon apithanos", "Pseudoloricaria laeviuscula", "Pseudorinelepis genibarbis", "Pterygoplichthys gibbiceps", "Rhinotocinclus eppleyi",
                         "Rhinelepis strigosa", "Rineloricaria lanceolata", "Scobinancistrus aureatus", "Spatuloricaria evansii", "Spectracanthicus murinus", "Spectracanthicus zuanoni", "Sturisoma robustum", "Sturisomatichthys tamanae")

treenooutgroup <- keep.tip(treeout, speciephynooutgroup)
dentaryfaceonnooutgroup <- readland.tps("dentary_front_faceon phy LM no outgroup.TPS", specID = "ID", negNA = TRUE, readcurves= TRUE)
dentfaceonnooutgroup <- gpagen(dentaryfaceonnooutgroup)
dentarysidenooutgroup <- readland.tps("dentary_side phy LM1 no outgroup.TPS", specID = "ID", negNA = TRUE, readcurves= TRUE)
dentsidenooutgroup <- gpagen(dentarysidenooutgroup)
              
####CVA
#side
sideoralcvs <-read.table("CVA ..., scores side oral jaw against classiifers.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

PlotsideoralCV1CV2 <- ggscatter(sideoralcvs, x = "CV1", y = "CV2", color = "Oraljaw", font.label = c("italic", "black"),
                                    palette = c("blue", "limegreen", "orange"), 
                                    ellipse = FALSE, ellipse.type = "convex", ellipse.alpha = 0.08, mean.point = FALSE, star.plot = TRUE, legend.title = "Oral jaw shape") +
  theme_classic() +
  labs(x ="CV1 (61.20% of variation)", y = "CV2 (38.80% of variation)") +
  theme(legend.position = "none") 
PlotsideoralCV1CV2

#faceon
faceonoralcvs <-read.table("CVA ..., scores faceon oral against classifier.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

PlotfaceonoralCV1CV2 <- ggscatter(faceonoralcvs, x = "CV1", y = "CV2", color = "Oral jaw", font.label = c("italic", "black"),
                                    palette = c("blue", "limegreen", "orange"), 
                                    ellipse = FALSE, ellipse.type = "convex", ellipse.alpha = 0.08, mean.point = FALSE, star.plot = TRUE, legend.title = "Oral jaw shape") +
  theme_classic() +
  labs(x ="CV1 (89.18% of variation)", y = "CV2 (10.82% of variation)") +
  theme(legend.position = "none") 
PlotfaceonoralCV1CV2

############################################################################
############################################################################
############################################################################
############################################################################

######pharyngeal jaws ###############
##top#################################################diet
pharyngealtopright <- readland.tps("top two pharyngeal jaws phy lm diet no outgroup.TPS", specID = "ID", negNA = TRUE, readcurves= TRUE)
pharytopright <- gpagen(pharyngealtopright)
Lori1 <- gpagen(pharyngealtopright)
plot(pharytopright)
phygdf <- geomorph.data.frame(pharytopright, categories = categories)
phyaov <- procD.pgls(phygdf$coords ~ categories, phy = tree) 
anova(phyaov)  #######           Df     SS      MS     Rsq     F      Z Pr(>F)  
                     #categories  2 179.92 89.96 0.79696 82.428 2.4301  0.01 **
                     #Residuals  42 45.838 1.091 0.20304                     
                     #Total      44 225.758                                      

toppharycvs <-read.table("CVA ..., scores for diet with two jaws.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

PlottoppharydietCVA <- ggscatter(toppharycvs, x = "CV1", y = "CV2", color = "Diet", font.label = c("italic", "black"),
                                 palette = c('Detritivore' = "blue", 'Algivore' = "limegreen", 'Insectivore' = "orange"),
                     ellipse = FALSE, ellipse.type = "convex", ellipse.alpha = 0.08, mean.point = FALSE, star.plot = TRUE, legend.title = "Dietary Niche") +
  theme_classic() +
  labs(x ="CV1 (74.25% of variation)", y = "CV2 (25.75% of variation)") +
  theme(legend.position = "right") 
PlottoppharydietCVA

######side############################################
#####pharyngeal jaw#########################################Diet
pharyngealside <- readland.tps("side_pharyngeal_jaw phylm diet no outgroup.TPS", specID = "ID", negNA = TRUE, readcurves= TRUE)
pharyside <- gpagen(pharyngealside)
Lori1 <- gpagen(pharyngealside)
plot(pharyside)
phygdf <- geomorph.data.frame(pharyside, categories = categories)
phyaov <- procD.pgls(phygdf$coords ~ categories, phy = tree) 
anova(phyaov)   ########            Df     SS      MS     Rsq      F      Z Pr(>F)  
                       #categories  2 37.002 18.5012 0.54007 24.659 2.5875  0.012 *
                       #Residuals  42 31.512  0.7503 0.45993                       
                       #Total      44 68.514                   

sidepharycvs <-read.table("CVA scores side pharyngeal jaw diet.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

PlotsidepharydietCVA <- ggscatter(sidepharycvs, x = "CV1", y = "CV2", color = "cat", font.label = c("italic", "black"),
                                  palette = c('Detritivore' = "blue", 'Algivore' = "limegreen", 'Insectivore' = "orange"),
                                  ellipse = FALSE, ellipse.type = "convex", ellipse.alpha = 0.08, mean.point = FALSE, star.plot = TRUE, legend.title = "Dietary Niche") +
  theme_classic() +
  labs(x ="CV1 (70.92% of variation)", y = "CV2 (29.08% of variation)") +
  theme(legend.position = "right") 

PlotsidepharydietCVA 


#####classifyier reliability
#######comparing pharyngeal jaw classifiers to pharyngeal jaw shape

toppharycvs <-read.table("CVA ..., scores top pharyngeal jaw to classifiers.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

PlottoppharyclassifyCVA <- ggscatter(toppharycvs, x = "CV1", y = "CV2", color = "pharyngeal jaws", font.label = c("italic", "black"),
                                  palette = c("blue", "limegreen", "orange", "purple", "yellow"),
                                  ellipse = FALSE, ellipse.type = "convex", ellipse.alpha = 0.08, mean.point = FALSE, star.plot = TRUE, legend.title = "Pharyngeal jaw shape") +
  theme_classic() +
  labs(x ="CV1 (70.92% of variation)", y = "CV2 (29.08% of variation)") +
  theme(legend.position = "right") 

PlottoppharyclassifyCVA 

##################################################
#######compared to oral jaws to identify link 
toppharynooutgroup <- readland.tps("top two pharyngeal jaws phy lm no outgroup.TPS", specID = "ID", negNA = TRUE, readcurves= TRUE)
toppharynooutgroup <- gpagen(toppharynooutgroup)

IDphary <- read.csv("IDpharyphy no outgroup.csv", header = T)
pharycat <-factor(IDphary$oral.jaw.1)

###top
phygdftopphary <- geomorph.data.frame(toppharynooutgroup, pharycat = pharycat)
phyaov <- procD.pgls(phygdftopphary$coords ~ pharycat, phy = treenooutgroup) 
anova(phyaov)   ##########           Df     SS     MS     Rsq      F     Z Pr(>F)   
                         #pharycat   2 178.520 89.260 0.72641 79.654 2.3477  0.001 ***
                         #Residuals 60  67.236  1.121 0.27359                         
                         #Total     62 245.756  



toppharyoralcvs <-read.table("CVA ..., scores from pharyngeal to oral with both jaws.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

PlotpharytoporalCV1CV2 <- ggscatter(toppharyoralcvs, x = "CV1", y = "CV2", color = "Oral jaw shape", font.label = c("italic", "black"),
                                     palette = c('Typical' ="black", 'reduced' ="magenta", 'elongate and long at toothcup' ="red"), 
                                     ellipse = FALSE, ellipse.type = "convex", ellipse.alpha = 0.08, mean.point = FALSE, star.plot = TRUE, legend.title = "Oral jaw shape") +
  theme_classic() +
  labs(x ="CV1 (78.53% of variation)", y = "CV2 (21.47% of variation)") 
 # theme(legend.position = "none") 
PlotpharytoporalCV1CV2



###combined figures

plotsjaws <- ggarrange(plotmds, dietcomponents, PlottoppharydietCVA, PlotpharytoporalCV1CV2,
                       nrow = 2, ncol =2,
                       labels = c("A.", "B.", "C.", "D.", "E.", "F."), align = "hv")
plotsjaws


