library(ggplot2)
library(plyr)
library(dplyr)
library("ggsci")
library("gridExtra")
library(grid)
library(RColorBrewer)
library("ggsci")
library("ggpubr")
library(stringr)

#display.brewer.all()
nb.cols <- 10
mycolors <- colorRampPalette(brewer.pal(10, "Set3"))(nb.cols)

##### Bacterial diversity plots

df <- read.table("../data/batch1_readcount_table_forR_noContaminants.tsv", sep="\t", h=T)
metadata <- read.csv("../data/metadata_macotra_batch1_wSAabund.tsv", sep="\t", h=T)

names(df) <- gsub(x = names(df), pattern = "\\.", replacement = "_")
metadata$sample <- gsub(x = metadata$sample, pattern = "\\-", replacement = "_")

Ncarrier <- c(as.character(metadata[metadata$carriage==0,]$sample))
carrier <- c(as.character(metadata[metadata$carriage==1,]$sample))

list_species <- c("Staphylococcus_aureus", "Staphylococcus_epidermidis", "Cutibacterium_acnes", "Corynebacterium_accolens", "Corynebacterium_macginleyi", "Corynebacterium_propinquum", "Corynebacterium_pseudodiphtheriticum", "Dolosigranulum_pigrum", "Moraxella_nonliquefaciens")

sub_species <- df[df$Species %in% list_species,]
sub_species = select(sub_species, -c("Genus","Family","Order","Class","Phylum"))
sub_others <- df[! df$Species %in% list_species,]
sub_others = select(sub_others, -c("Species","Genus","Family","Order","Class","Phylum"))
sub_others_sum <- colSums(sub_others)
sub_others_sum <- as.data.frame(t(sub_others_sum))
sub_others_sum$Species <- "Others"
df <- rbind(sub_species, sub_others_sum)


Ncarrier_data <- metadata[metadata$carriage==0,]
carrier_data <- metadata[metadata$carriage==1,]

data_Ncarrier <- df[,c(Ncarrier,"Species")]
data_carrier <- df[,c(carrier,"Species")]



#Carrier
T0_id <- as.vector(carrier_data[carrier_data$sample_time=="T0",]$sample)
T1_id <- as.vector(carrier_data[carrier_data$sample_time=="T1",]$sample)
T2_id <- as.vector(carrier_data[carrier_data$sample_time=="T2",]$sample)
T3_id <- as.vector(carrier_data[carrier_data$sample_time=="T3",]$sample)
T4_id <- as.vector(carrier_data[carrier_data$sample_time=="T4",]$sample)

T0_data <- data_carrier[,T0_id]
T1_data <- data_carrier[,T1_id]
T2_data <- data_carrier[,T2_id]
T3_data <- data_carrier[,T3_id]
T4_data <- data_carrier[,T4_id]

T0_data_to_use <- data.frame(Species=data_carrier$Species, Means=rowMeans(T0_data))
T1_data_to_use <- data.frame(Species=data_carrier$Species, Means=rowMeans(T1_data))
T2_data_to_use <- data.frame(Species=data_carrier$Species, Means=rowMeans(T2_data))
T3_data_to_use <- data.frame(Species=data_carrier$Species, Means=rowMeans(T3_data))
T4_data_to_use <- data.frame(Species=data_carrier$Species, Means=rowMeans(T4_data))

T0_data_to_use$time <- "J0"
T1_data_to_use$time <- "J7"
T2_data_to_use$time <- "M1"
T3_data_to_use$time <- "M3"
T4_data_to_use$time <- "M6"

df_graph <- rbind(T0_data_to_use, T1_data_to_use, T2_data_to_use, T3_data_to_use, T4_data_to_use)

df_graph$Species <- factor(df_graph$Species, levels = c("Staphylococcus_aureus", "Dolosigranulum_pigrum", "Moraxella_nonliquefaciens", "Corynebacterium_propinquum", "Corynebacterium_accolens", "Corynebacterium_pseudodiphtheriticum", "Corynebacterium_macginleyi", "Staphylococcus_epidermidis", "Cutibacterium_acnes", "Others")
)

p1 <- ggplot(data=df_graph, aes(x=time, y=Means, fill=Species)) +
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values=mycolors) + theme_minimal() + ylab("Proportions") + xlab("Sampling time") +
  scale_x_discrete(limits=c("J0", "J7", "M1", "M3", "M6")) +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.title = element_text(size = 18),
        legend.position = "none")
p1 <- p1 + geom_vline(xintercept = 1.5, color="red",linetype = 2, size = 0.7) 
#p1 <- p1 + theme(axis.title = element_text(size= 20 ))
p1



#Ncarrier
T0_id <- as.vector(Ncarrier_data[Ncarrier_data$sample_time=="T0",]$sample)
T1_id <- as.vector(Ncarrier_data[Ncarrier_data$sample_time=="T1",]$sample)
T2_id <- as.vector(Ncarrier_data[Ncarrier_data$sample_time=="T2",]$sample)
T3_id <- as.vector(Ncarrier_data[Ncarrier_data$sample_time=="T3",]$sample)
T4_id <- as.vector(Ncarrier_data[Ncarrier_data$sample_time=="T4",]$sample)

T0_data <- data_Ncarrier[,T0_id]
T1_data <- data_Ncarrier[,T1_id]
T2_data <- data_Ncarrier[,T2_id]
T3_data <- data_Ncarrier[,T3_id]
T4_data <- data_Ncarrier[,T4_id]

T0_data_to_use <- data.frame(Species=data_Ncarrier$Species, Means=rowMeans(T0_data))
T1_data_to_use <- data.frame(Species=data_Ncarrier$Species, Means=rowMeans(T1_data))
T2_data_to_use <- data.frame(Species=data_Ncarrier$Species, Means=rowMeans(T2_data))
T3_data_to_use <- data.frame(Species=data_Ncarrier$Species, Means=rowMeans(T3_data))
T4_data_to_use <- data.frame(Species=data_Ncarrier$Species, Means=rowMeans(T4_data))

T0_data_to_use$time <- "J0"
T1_data_to_use$time <- "J7"
T2_data_to_use$time <- "M1"
T3_data_to_use$time <- "M3"
T4_data_to_use$time <- "M6"

df_graph <- rbind(T0_data_to_use, T1_data_to_use, T2_data_to_use, T3_data_to_use, T4_data_to_use)


#df_graph$Species <- factor(df_graph$Species, levels = c("Staphylococcus_aureus", "Dolosigranulum_pigrum", "Moraxella_nonliquefaciens", "Corynebacterium_propinquum", "Corynebacterium_accolens", "Corynebacterium_pseudodiphtheriticum", "Corynebacterium_macginleyi", "Staphylococcus_epidermidis", "Cutibacterium_acnes", "Others"))
selected_names <- c("Staphylococcus_aureus", "Dolosigranulum_pigrum", "Moraxella_nonliquefaciens", "Corynebacterium_propinquum", "Corynebacterium_accolens", "Corynebacterium_pseudodiphtheriticum", "Corynebacterium_macginleyi", "Staphylococcus_epidermidis", "Cutibacterium_acnes", "Others")
df_graph$Species <- factor(df_graph$Species, levels = selected_names)
clean_names <- c("S. aureus", "D. pigrum", "M. nonliquefaciens", "C. propinquum", "C. accolens", "C. pseudodiphtheriticum", "C. macginleyi", "S. epidermidis", "C. acnes", "Others")
df_graph$Species <- df_graph$Species %>% str_replace_all(c("Staphylococcus_aureus" = "S. aureus", 
                                                           "Dolosigranulum_pigrum" = "D. pigrum", 
                                                           "Moraxella_nonliquefaciens" = "M. nonliquefaciens",
                                                           "Corynebacterium_propinquum" = "C. propinquum", 
                                                           "Corynebacterium_accolens" = "C. accolens", 
                                                           "Corynebacterium_pseudodiphtheriticum" = "C. pseudodiphtheriticum", 
                                                           "Corynebacterium_macginleyi" = "C. macginleyi", 
                                                           "Staphylococcus_epidermidis" = "S. epidermidis", 
                                                           "Cutibacterium_acnes" = "C. acnes"))
df_graph$Species <- factor(df_graph$Species, levels = clean_names)

# p2 <- ggplot(data=df_graph, aes(x=time, y=Means, fill=Species)) +
#   geom_bar(stat="identity", position="fill") +
#   scale_fill_manual(values=mycolors) + theme_minimal() + ylab("Proportions") + xlab("Days after decolonization") + 
#   ggtitle("Bacterial abundance in non-carriers") + scale_x_discrete(limits=c("0", "7", "35", "91", "175")) + 
#   theme(legend.text = element_text(size=18), 
#         legend.title = element_text(size=18),
#         axis.text.x = element_text(size = 15),
#         axis.text.y = element_text(size = 15),
#         axis.title.x = element_text(size = 15),
#         axis.title.y = element_text(size = 15),
#         plot.title = element_text(hjust = 0.5, face= "bold")) 
# p2 <- p2 + geom_vline(xintercept = 1.5, color="red",linetype = 2, size = 0.5)
# p2

FUNitalize<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Others", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

p2 <- ggplot(data=df_graph) +
  geom_bar(aes(x=time, y=Means, fill=Species), stat="identity", position="fill") +
  scale_fill_manual(values=mycolors, breaks=levels(df_graph$Species), labels = FUNitalize(levels(df_graph$Species))) + theme_minimal() + ylab("Proportions") + xlab("Sampling time") + 
  scale_x_discrete(limits=c("J0", "J7", "M1", "M3", "M6")) + 
  theme(legend.text = element_text(size=17), 
        legend.text.align = 0,
        legend.title = element_text(size=20),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.title = element_text(size = 18))

p2 <- p2 + geom_vline(xintercept = 1.5, color="red",linetype = 2, size = 0.7)
p2


grid.arrange(p1,p2, widths = c(4,10) )

svg("../outputs/figure_bacterial_abund.svg", width=1200, height=796)
# Note : Zoom -> save image
dev.off()



##### Alpha diversity boxplots

df <- read.table("d_diversity.tsv", sep="\t", h=T)
df$sample <- gsub(x = df$sample, pattern = "\\-", replacement = "_")

data_Ncarrier <- filter(df, sample %in% Ncarrier)
data_carrier <- filter(df, sample %in% carrier)


T0_id <- as.vector(as.character(Ncarrier_data[Ncarrier_data$sample_time=="T0",]$sample))
T1_id <- as.vector(Ncarrier_data[Ncarrier_data$sample_time=="T1",]$sample)
T2_id <- as.vector(Ncarrier_data[Ncarrier_data$sample_time=="T2",]$sample)
T3_id <- as.vector(Ncarrier_data[Ncarrier_data$sample_time=="T3",]$sample)
T4_id <- as.vector(Ncarrier_data[Ncarrier_data$sample_time=="T4",]$sample)

T0_data <- filter(data_Ncarrier, sample %in% T0_id)$Species.shannon.Estimator
T1_data <- filter(data_Ncarrier, sample %in% T1_id)$Species.shannon.Estimator
T2_data <- filter(data_Ncarrier, sample %in% T2_id)$Species.shannon.Estimator
T3_data <- filter(data_Ncarrier, sample %in% T3_id)$Species.shannon.Estimator
T4_data <- filter(data_Ncarrier, sample %in% T4_id)$Species.shannon.Estimator


ncarrier_shannon <- data.frame(shannon_index=c(T0_data, T1_data, T2_data, T3_data, T4_data), Time=c(rep("T0", length(T0_data)), rep("T1", length(T1_data)), rep("T2", length(T2_data)), rep("T3", length(T3_data)), rep("T4", length(T4_data))))

my_comparaison <- list(c("T0", "T1"),  c("T0", "T2"), c("T0", "T3"), c("T0", "T4"), c("T1", "T2"), c("T1", "T3"), c("T1", "T4"), c("T2", "T3"), c("T3", "T4"))
p3 <- ggboxplot(ncarrier_shannon, x= "Time", y="shannon_index", color = "Time") + stat_compare_means(comparisons = my_comparaison, aes(label=..p.adj..)) + ggtitle("Non carriers")



T0_id <- as.vector(carrier_data[carrier_data$sample_time=="T0",]$sample)
T1_id <- as.vector(carrier_data[carrier_data$sample_time=="T1",]$sample)
T2_id <- as.vector(carrier_data[carrier_data$sample_time=="T2",]$sample)
T3_id <- as.vector(carrier_data[carrier_data$sample_time=="T3",]$sample)
T4_id <- as.vector(carrier_data[carrier_data$sample_time=="T4",]$sample)

T0_data <- filter(data_carrier, sample %in% T0_id)$Species.shannon.Estimator
T1_data <- filter(data_carrier, sample %in% T1_id)$Species.shannon.Estimator
T2_data <- filter(data_carrier, sample %in% T2_id)$Species.shannon.Estimator
T3_data <- filter(data_carrier, sample %in% T3_id)$Species.shannon.Estimator
T4_data <- filter(data_carrier, sample %in% T4_id)$Species.shannon.Estimator


carrier_shannon <- data.frame(shannon_index=c(T0_data, T1_data, T2_data, T3_data, T4_data), Time=c(rep("T0", length(T0_data)), rep("T1", length(T1_data)), rep("T2", length(T2_data)), rep("T3", length(T3_data)), rep("T4", length(T4_data))))

my_comparaison <- list(c("T0", "T1"),  c("T0", "T2"), c("T0", "T3"), c("T0", "T4"), c("T1", "T2"), c("T1", "T3"), c("T1", "T4"), c("T2", "T3"), c("T3", "T4"))
p4 <- ggboxplot(carrier_shannon, x= "Time", y="shannon_index", color = "Time") + stat_compare_means(comparisons = my_comparaison, aes(label=..p.adj..)) + ggtitle("Carriers")
p4




#pdf("figure_article_shannon.pdf", width=30, height=18)
#grid.arrange(p2,p4,p1,p3, ncol=2, nrow=2,top = textGrob("Carriers vs Non carriers"))


