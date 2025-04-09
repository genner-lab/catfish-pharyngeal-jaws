#############Pharyngeal or excluding most of Ancistrus Jaw script (For all Loricariidae)#########
library(geomorph) #
library(ape) #
library(ggplot2) #
library(MASS) # 
library(ggpubr) #
setwd("C:/Users/sl22263/OneDrive - University of Bristol/Documents/Data to be landmarked/CT scan images/Loricariidae Samples/analyses")
#########jaw type identifyers
ID <- read.csv("jawtype.csv",header=TRUE,fill=TRUE,sep=",",check.names=FALSE)
Species <- factor(ID$Corrected)
phary <- factor(ID$pharyngealjaw)
oralinfo1 <- factor(ID$`oral jaw 1`)

#####oral jaws#########
dentaryfaceon <- readland.tps("dentary_front_faceonphy.TPS", specID = "ID", negNA = TRUE, readcurves= TRUE)
dentfaceon <- gpagen(dentaryfaceon)
dentaryside <- readland.tps("dentary_sidephy.TPS", specID = "ID", negNA = TRUE, readcurves= TRUE)
dentside <- gpagen(dentaryside)
dentaryfrontbelow <- readland.tps("dentary_front_belowphy.TPS", specID = "ID", negNA = TRUE, readcurves= TRUE)
dentfrontbelow <- gpagen(dentaryfrontbelow)
dentarytop <- readland.tps("dentarytopphy.TPS", specID = "ID", negNA = TRUE, readcurves= TRUE)
denttop <- gpagen(dentarytop)
maxilliae <- readland.tps("maxilliae_fronttest.TPS", specID = "ID", negNA = TRUE, readcurves= TRUE)
maxil <- gpagen(maxilliae)
oral <- combine.subsets(faceon = dentfaceon, side = dentside, below = dentfrontbelow, top = denttop, max = maxil,  gpa = TRUE)


summary(oral)
plotAllSpecimens(oral$coords)
oralshape <- two.d.array(oral$coords)
oraldata <- as.data.frame(oralshape)
#pca
PCAoral <- prcomp(oraldata[2:1000], scale = TRUE)
Example_PCAoral_score <- as.data.frame(PCAoral$x)
PCAoral_scores <- cbind(oralinfo,Example_PCAoral_score)
#lda
Plecooral_LDA_Data <- cbind(ID,oraldata[2:1000])
LDAoral <- lda(Plecooral_LDA_Data[,c(7,1004)], Plecooral_LDA_Data$`oral jaw 1`)
LDAoral_predict <- predict(LDAoral)
LDAoral_score <- as.data.frame(LDAoral_predict$x)
LDAoral_score <- cbind(ID[1:5], LDAoral_score)
PlotLDA <- ggscatter(LDAoral_score, x = "LD1", y = "LD2", color = "oral jaw 1", label = "Corrected", font.label = c("italic", "black"),
                     ellipse = FALSE, ellipse.type = "convex", ellipse.alpha = 0.08, mean.point = FALSE, star.plot = TRUE, legend.title = "Dietary Niche") + 
  theme_classic() +
  labs(x ="DF1 (X% of variation)", y = "DF2 (X% of variation)") +
  theme(legend.position = "right") 
PlotLDA


######pharyngeal jaws ###############
pharyngealtop <- readland.tps("top_pharyngeal_jawphy.TPS", specID = "ID", negNA = TRUE, readcurves= TRUE)
pharytop <- gpagen(pharyngealtop)
pharyngealside <- readland.tps("side_pharyngeal_jawphy.TPS", specID = "ID", negNA = TRUE, readcurves= TRUE)
pharyside <- gpagen(pharyngealside)
pharyngeal <- combine.subsets(side = pharyside, top = pharytop, gpa = TRUE)
summary(pharyngeal)
plotAllSpecimens(pharyngeal$coords)
pharyngealshape <- two.d.array(pharyngeal$coords)
pharyngealdata <- as.data.frame(pharyngealshape)
#pca
PCApharyngeal <- prcomp(pharyngealdata[2:400], scale = TRUE)
Example_PCAphary_score <- as.data.frame(PCApharyngeal$x)
PCAphary_scores <- cbind(phary,Example_PCAphary_score)
PlotPCA <- ggscatter(PCAphary_scores, x = "PC1", y = "PC2", color = "phary", font.label = c("italic", "black"),
                     ellipse = FALSE, ellipse.type = "convex", ellipse.alpha = 0.08, mean.point = FALSE, star.plot = TRUE, legend.title = "Dietary Niche") + 
  theme_classic() +
  geom_point(aes(colour= "Molariform")) +
  labs(x ="PC1 (X% of variation)", y = "PC2 (X% of variation)") +
  theme(legend.position = "right") 
PlotPCA

########GLM 
tree <- read.nexus(file = "phylogenyLoricariidae2.nex")
library(phylolm)
ID <- read.csv("jawtypeglm.csv",header=TRUE,fill=TRUE,sep=",",check.names=FALSE)
Jaw <- factor(ID$test)
jaw = data.frame(trait01 = ID$test, predictor = ID$`oral jaw 1`)
rownames(jaw) <- ID$Corrected
fit = phyloglm(trait01~predictor, phy=tree,data=jaw,boot=1000)
summary(fit)
coef(fit)
##alpha: 0.07164009 
#bootstrap mean: 0.4418443 (on log scale, then back transformed)
#so possible upward bias.
#bootstrap 95% CI: (0.04423241,41.53958)
#Coefficients:
 # Estimate      StdErr     z.value lowerbootCI upperbootCI p.value
#(Intercept)                            -1.1519e+00  2.3394e+00 -4.9239e-01 -2.2687e+00      0.1996  0.6224
#predictorelongate and long at toothcup -1.9044e-02  4.8611e-01 -3.9176e-02 -1.6505e+00      2.2016  0.9687
#predictorreduced                        2.8961e+00  1.9607e+00  1.4771e+00 -4.4288e-14      4.7449  0.1397
#predictorTypical                       -6.9080e-01  6.0961e-01 -1.1332e+00 -3.3146e+00      0.0551  0.2571
###########do I need btol?

#########Contingency table ##################sort out###########
jawdf <- data.frame(ID$test, ID$`oral jaw 1`)
jawtab <- table(jawdf)
jawtab
tab <- mosaicplot(jawtab)

##diet with pharyngeal
#########classifiers
ID <- read.csv("IDpharyphy.csv", header = T)
Genus <- factor(ID$Genus)
Diet <- factor(ID$Niche)
Subfamily <-factor(ID$Subfamily)
Location <-factor(ID$Substrate)
Species <-factor(ID$name)
##########ANOVA phary side
phygdf <- geomorph.data.frame(pharyside, diet = Diet)
phyaov <- procD.pgls(phygdf$coords ~ Diet, phy = tree) 
anova(phyaov)
##Estimation method: Generalized Least-Squares (via OLS projection) 
####Sums of Squares and Cross-products: Type I 
#####Effect sizes (Z) based on F distributions
#########ANOVA phary top
plot(phyaov)
phygdf <- geomorph.data.frame(pharytop, diet = Diet)
phyaov <- procD.pgls(phygdf$coords ~ Diet, phy = tree) 
anova(phyaov)
plot(phyaov)
##Estimation method: Generalized Least-Squares (via OLS projection) 
####Sums of Squares and Cross-products: Type I 
#####Effect sizes (Z) based on F distributions

