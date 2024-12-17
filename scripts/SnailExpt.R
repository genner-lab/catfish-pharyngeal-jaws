setwd("C:/Users/sl22263/OneDrive - University of Bristol/Documents/Data to be landmarked/CT scan images/Loricariidae Samples/analyses/Pharyngeal figures/Feeding trial")
Snail_Expt <- read.table("Snails3.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)

library(lme4)
library(lmerTest)
library(emmeans)
library (ggplot2)
library(ggpubr)

Model1 <- lmer(Small_consumed ~ Species + (1|FishID) + (1|TrialDay), data = Snail_Expt)
summary(Model1)
Model1_means <- emmeans(Model1, "Species")

Model2 <- lmer(Large_consumed ~ Species + (1|FishID) + (1|TrialDay), data = Snail_Expt)
summary(Model2)
Model2_means <- emmeans(Model2, "Species")

Model3 <- lmer(Malaysian_consumed ~ Species + (1|FishID) + (1|TrialDay), data = Snail_Expt)
summary(Model3)
Model3_means <- emmeans(Model3, "Species")

Snail_Data2 <- read.table("Snails3.txt",header=TRUE,fill=TRUE,sep="\t",check.names=FALSE)



p1<- ggplot(Snail_Data2, aes(x=FishID, y=Small_consumed, fill=Species)) + 
  geom_bar(stat="summary", color="black", 
           position=position_dodge()) + theme_classic() + ylim(0,10) +
  scale_fill_manual(values=c("#DDAA33", "#4477AA")) +
  labs(x = "Individual fish", y = "Mean number of small pulmonate consumed")
p1

p2<- ggplot(Snail_Data2, aes(x=FishID, y=Large_consumed, fill=Species)) + 
  geom_bar(stat="summary", color="black", 
           position=position_dodge()) + theme_classic() + ylim(0,10) +
  scale_fill_manual(values=c("#DDAA33", "#4477AA")) +
  labs(x = "Individual fish", y = "Mean number of large pulmonate consumed")
p2
p3<- ggplot(Snail_Data2, aes(x=FishID, y=Malaysian_consumed, fill=Species)) + 
  geom_bar(stat="summary", color="black", 
           position=position_dodge()) + theme_classic() + ylim(0,10) +
  scale_fill_manual(values=c("#DDAA33", "#4477AA")) +
  labs(x = "Individual fish", y = "Mean number of prosobranch consumed")
p3

p_both <- ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = TRUE, legend = "right")
p_both


