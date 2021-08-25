library(data.table)
library(readxl)
library(iNEXT)
library(dplyr)
library(tidyverse)
library("gridExtra")
library(vegan)


# Distance between population structures

d <- fread("../data/batch1_readcount_table_forR_noContaminants.tsv", sep ="\t", h=T)
m <- fread("../data/metadata_macotra_batch1_wSAabund.tsv", sep ="\t", h=T)

list_species <- c("Staphylococcus_aureus", "Staphylococcus_epidermidis", "Cutibacterium_acnes", "Corynebacterium_accolens", "Corynebacterium_macginleyi", "Corynebacterium_propinquum", "Corynebacterium_pseudodiphtheriticum", "Dolosigranulum_pigrum", "Moraxella_nonliquefaciens")

sub_species <- d[d$Species %in% list_species,]
sub_species = select(sub_species, -c("Genus","Family","Order","Class","Phylum"))
sub_others <- d[! d$Species %in% list_species,]
sub_others = select(sub_others, -c("Species","Genus","Family","Order","Class","Phylum"))
sub_others_sum <- colSums(sub_others)
sub_others_sum <- as.data.frame(t(sub_others_sum))
sub_others_sum$Species <- "Others"
df <- rbind(sub_species, sub_others_sum)
df2 <- as.data.table(head(t(df), -1))
colnames(df2) <- c(list_species, "Others")
df3 <- df2 %>% mutate_if(is.character, ~as.numeric(as.character(.)))
all(unlist(lapply(df3, class)) == "numeric")
#df3 <- df3 / rowSums(df3)
df4 <- cbind(m[, .(sample, carriage, participant_id, days_after_start_of_treatment)], df3)



bray <- vegdist(df3)

hist(bray)
range(bray)

Ncarrier_list <- levels(factor(df4[df4$carriage==0,]$participant_id))
carrier_list <- levels(factor(df4[df4$carriage==1,]$participant_id))

days <- levels(factor(df3$days_after_start_of_treatment))
patients <- levels(factor(df3$participant_id))

braymat <- as.matrix(bray)
jind <- rep(1 + 5*(0:15), each = 5)
braymat2 <- braymat[1:nrow(braymat), jind]
braydiv <- df4 %>% select(days_after_start_of_treatment, participant_id)
braydiv$bray_curtis_toT0 <- diag(braymat2)
braydiv

#mean_BC <- group_by(braydiv2, days_after_start_of_treatment) %>% summarise(mean_BC_toT0 = mean(bray_curtis_toT0, na.rm = TRUE))

braydiv_carrier <- subset(braydiv, participant_id %in% carrier_list)
braydiv_Ncarrier <- subset(braydiv, participant_id %in% Ncarrier_list)

braydiv_carrier_mean <- braydiv_carrier[, .(BC_mean = mean(bray_curtis_toT0), BC_cilo = mean(bray_curtis_toT0) - 1.96*sd(bray_curtis_toT0)/sqrt(8), BC_cihi = mean(bray_curtis_toT0) + 1.96*sd(bray_curtis_toT0)/sqrt(8)), by = .(days_after_start_of_treatment)  ]
braydiv_Ncarrier_mean <- braydiv_Ncarrier[, .(BC_mean = mean(bray_curtis_toT0), BC_cilo = mean(bray_curtis_toT0) - 1.96*sd(bray_curtis_toT0)/sqrt(9), BC_cihi = mean(bray_curtis_toT0) + 1.96*sd(bray_curtis_toT0)/sqrt(9)), by = .(days_after_start_of_treatment)  ]

p1 <- ggplot(braydiv_carrier_mean, aes(x = days_after_start_of_treatment, y = BC_mean)) + geom_point(colour = 'red', size = 4, shape = 4) + geom_line(colour = 'red', size = 1) +
  ylab("Bray-Curtis dissimilarity to day 0") + xlab("Days after decolonization") + ylim(0,1) + xlim(0, 175) + 
  ggtitle("Dynamics of Bray-Curtis dissimilarities compared to day 0 in carriers") + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_ribbon(aes(ymin=BC_cilo, ymax=BC_cihi), linetype=2, alpha=0.1)

for(participant_current in levels(factor(braydiv_carrier$participant_id))) {
  dsub <- braydiv_carrier[participant_id == participant_current]
  p1 <- p1 + geom_line(data = dsub,aes(days_after_start_of_treatment,bray_curtis_toT0), linetype = 3)
}
p1


p2 <- ggplot(braydiv_Ncarrier_mean, aes(x = days_after_start_of_treatment, y = BC_mean)) + geom_point(colour = 'red', size = 4, shape = 4) + geom_line(colour = 'red', size = 1) +
  ylab("Bray-Curtis dissimilarity to day 0") + xlab("Days after decolonization") + ylim(0,1) + xlim(0, 175) + 
  ggtitle("Dynamics of Bray-Curtis dissimilarities compared to day 0 in non carriers") + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_ribbon(aes(ymin=BC_cilo, ymax=BC_cihi), linetype=2, alpha=0.1)

for(participant_current in levels(factor(braydiv_Ncarrier$participant_id))) {
  dsub <- braydiv_Ncarrier[participant_id == participant_current]
  p2 <- p2 + geom_line(data = dsub,aes(days_after_start_of_treatment,bray_curtis_toT0), linetype = 3)
}
p2

svg("../outputs/Bray-curtis_dynamics.svg", width=1300, height=500)
grid.arrange(p1,p2, ncol=2, nrow=1)
#Note : manually export to svg because bug.
dev.off()
