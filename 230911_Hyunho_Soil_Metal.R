# Alpha diversity
library(ggplot2)
library(reshape2)
library(dplyr)
# Font
library(extrafont)
font_import()
fonts()
loadfonts(device="win")

setwd("D:/220729_Metagenomics_R_script/220808_R/230911_Hyunho_Soil_Metal/230912_Alpha_diversity/Alpha_data")
getwd()

A <- read.table("230912_Alpha_diversity.txt", sep = "\t", header = TRUE)
A$Treatment <- gsub("observed_features","Observed ASVs",A$Treatment)
A$Treatment <- gsub("faith_pd","Faith's phylogenetic diversity",A$Treatment)
A$Treatment <- gsub("simpson","Simpson index",A$Treatment)
A$Treatment <- gsub("pielou_evenness","Evenness",A$Treatment)
A$Treatment <- gsub("shannon_entropy","Shannon index",A$Treatment)
A$Treatment <- gsub("chao1","Chao1",A$Treatment)
A$Treatment <- factor(A$Treatment,levels=c("Observed ASVs","Chao1","Evenness","Faith's phylogenetic diversity",
                                           "Shannon index","Simpson index"))
A$Treatment
A$variable <- factor(A$variable, levels = c("CON","MS","ZVM"))

B <- read.table("230912_Alpha_P_value.txt", sep = "\t", header = TRUE)
B$Treatment <- gsub("observed_features","Observed ASVs",B$Treatment)
B$Treatment <- gsub("faith_pd","Faith's phylogenetic diversity",B$Treatment)
B$Treatment <- gsub("simpson","Simpson index",B$Treatment)
B$Treatment <- gsub("pielou_evenness","Evenness",B$Treatment)
B$Treatment <- gsub("shannon_entropy","Shannon index",B$Treatment)
B$Treatment <- gsub("chao1","Chao1",B$Treatment)
B$Treatment <- factor(B$Treatment,levels=c("Observed ASVs","Chao1","Evenness","Faith's phylogenetic diversity",
                                           "Shannon index","Simpson index"))
B

TRT_color <- c("CON"="#f6d55c", "MS"="#FFA080","ZVM"="#80A0FF")
# font --> face = 1 ---
par(mfrow=c(2,2))
name=c("1 : plane","2 : bold","3 : italic","4 : bold and italic")

# Box plot
g = ggplot(A, aes(x = variable, y = value, fill = variable, color=variable)) + 
  geom_boxplot(colour = "black", position = position_dodge(0.5), alpha=0.6, outlier.shape = NA) +
  stat_summary(fun=mean, geom="point", shape=17, size=3, color="black", fill="black") +
  geom_vline(xintercept = c(1.5,2.5,3.5), colour = "white", size = 0.8) +
  #geom_jitter(size=1.5, alpha=0.8)+
  theme(legend.title = element_blank(), 
        legend.text = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(face = "bold",colour = "black", size = 20, family = "Times new roman"), 
        axis.text.y = element_text(face = "bold", size = 20, colour = "black", family = "Times new roman"), 
        axis.title.y = element_text(face = "bold", size = 22, colour = "black", family = "Times new roman"),
        strip.text.x = element_text(face = "bold", size = 20, colour = "black", family = "Times new roman"))+ 
  labs(x= "", y = "Alpha diversity", fill = "") + 
  scale_fill_manual(values =  TRT_color, label = c("CON","MS","ZVM"))+
  scale_color_manual(values =  TRT_color, label = c("CON","MS","ZVM"))
g
g1 <- g+ facet_wrap(. ~ Treatment, scales = "free", ncol=2) +
  geom_text(data = B, aes(x=3.1, y=y, label = Label),
            size = 6,family = "Times New Roman", fontface="bold",
            colour = "black", inherit.aes=FALSE, parse = TRUE)
g1
ggsave("230912_Alpha_diversity_Box_plot.png", width = 9.5, height = 11, dpi = 600)


# Violin plot
g = ggplot(A, aes(x = variable, y = value, fill = variable, color=variable)) + 
  geom_violin(colour = "black", position = position_dodge(0.5)) +
  geom_vline(xintercept = c(1.5,2.5,3.5), colour = "white", size = 0.8) +
  #geom_jitter(size=0.8, alpha=0.9)+
  theme(legend.title = element_blank(), 
        legend.text = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(face = "bold",colour = "black", size = 20, family = "Times new roman"), 
        axis.text.y = element_text(face = "bold", size = 20, colour = "black", family = "Times new roman"), 
        axis.title.y = element_text(face = "bold", size = 22, colour = "black", family = "Times new roman"),
        strip.text.x = element_text(face = "bold", size = 20, colour = "black", family = "Times new roman"))+ 
  labs(x= "", y = "Alpha diversity", fill = "") + 
  scale_fill_manual(values =  TRT_color, label = c("HF-Low","HF-High","LF-High"))+
  scale_color_manual(values =  TRT_color, label = c("HF-Low","HF-High","LF-High"))
g1 <- g+ facet_wrap(. ~ Treatment, scales = "free", ncol=2)
g1

# https://r-charts.com/distribution/violin-plot-mean-ggplot2/
## How to add average on the violin plot
### geom= "point", "crossbar","
g2 <- g1+ stat_summary(fun="mean",geom="point",color="black")
g2
g3 <- g1+ stat_summary(fun.data="mean_cl_boot",geom="pointrange",color="black")
g3

ggsave("221023_Alpha_diversity_violin_plot_Silva.png", width = 10, height = 10, dpi = 300)

# Alpha rarefaction curve
library('dplyr')
library(plyr)
library(scales)
library(reshape2)
getwd()
setwd("D:/220729_Metagenomics_R_script/220808_R/230911_Hyunho_Soil_Metal/230912_Alpha_diversity/Alpha_rarefaction/")

A <- read.table('230912_Alpha_rarefaction_Input.txt',sep = "\t",header = TRUE)
View(A)
AM <- melt(A, id.vars = c("Index", "Depth", "Iter"))
View(AM)
df1 <- ddply(AM, c("Index", "Depth", "variable"), summarise, mean=mean(value), sd=sd(value), n= length(value), se=sd/sqrt(n))
View(df1)
df1$Treatment <- gsub("\\_[1-3]", "", df1$variable)
View(df1)
colnames(df1)[4]="value"
df2 <- ddply(df1, c("Index", "Depth", "Treatment"), summarise, mean=mean(value), sd=sd(value), n= length(value), se=sd/sqrt(n))
unique(df2$Index)
colnames(df2)[4]="Value"
df2
df2$Index <- gsub("Faith_PD", "Faith's phylogenetic diversity", df2$Index)
df2$Index <- factor(df2$Index,levels=c("Observed ASVs","Chao1","Evenness","Faith's phylogenetic diversity",
                                     "Shannon index","Simpson index"))
df2$Index
df2$Treatment <- factor(df2$Treatment, levels = c("CON","MS","ZVM"))
df2$Treatment
# geom_line --> color='color', lwd='line size' linetype='twodash','solid','longdash'
# stat_smooth(method='lm') # https://kuduz.tistory.com/1118
B = ggplot(data = df2, aes(Depth, Value, color=Treatment, shape=Treatment))+
  geom_point(lwd=2.5)+
  geom_line(lwd=1.2, linetype='solid', stat = "identity")+
  #coord_cartesian(ylim = c(0,1700))+
  scale_x_continuous(breaks=seq(0,31120,5000))+
  #scale_y_continuous(breaks=seq(0,1700,100))+
  geom_errorbar(aes(ymin=Value-se, ymax=Value+se),
                position=position_dodge(0.5), width=800, lwd=0.2, color = 'black')+
  scale_color_manual(name="",
                     limits=c("CON", "MS", "ZVM"), # order of the label
                     labels = c(expression(paste(bold("CON"))), # relabel the label
                                expression(paste(bold("MS"))),
                                expression(paste(bold("ZVM")))),
                     values = c("#f6d55c", "#FFA080", "#80A0FF"))+ # label color
  scale_shape_manual(name="", limits=c("CON", "MS", "ZVM"), # order of the label
                     labels = c(expression(paste(bold("CON"))), # relabel the label
                                expression(paste(bold("MS"))),
                                expression(paste(bold("ZVM")))),
                     values = c(16,16,16))+
  labs(x = "Sequencing depth", y = "", family = "Times New Roman", size = "15", face = "bold")+
  theme(legend.position = "bottom")
B
B1 <- B+ facet_wrap(. ~ Index, scales = "free", ncol=2) +
  theme(axis.title = element_text(size=20, family = "Times New Roman", face="bold", color="black"),
        axis.title.y = element_text(size = 20, family = "Times New Roman", face = "bold"), 
        legend.title = element_text(size = 20, family = "Times New Roman", face = "bold", colour = "black", hjust = 0.5), 
        legend.text = element_text(size = 20, family = "Times New Roman", face = "bold", colour = "black", hjust = 0), 
        axis.text.y = element_text(colour = "black", family = "Times New Roman", size = 20, face = "bold"), 
        axis.text.x = element_text(colour = "black", family = "Times New Roman", size = 20, face = "bold"),
        strip.text.x = element_text(face = "bold", size = 20, colour = "black", family = "Times new roman"))
B1
ggsave("230912_Alpha_rarefaction.png", width = 12, height = 12, dpi = 300)


# Beta diversity
# https://github.com/jbisanz/qiime2R
## Install qiime2R
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")
devtools::install_github("gdmcdonald/notly")
install.packages("tidyverse")
# https://plotly.com/r/static-image-export/
install.packages("plotly")
library(plotly)
library(tidyverse)
library(qiime2R)
library(dplyr)
library(vegan)
library(gridExtra)
library(notly)

# Set working directory
getwd()
setwd("D:/220729_Metagenomics_R_script/220808_R/230911_Hyunho_Soil_Metal/230912_Beta_diversity/")

# Data importing from qza result file of beta diversity (pcoa_results.qza)
metadata<-read_q2metadata("../230911_Meta_data.txt")
Unwe_uni_frac<-read_qza("unweighted_unifrac_pcoa_results.qza")
we_uni_frac<-read_qza("weighted_unifrac_pcoa_results.qza")

# 2D Plot
color <- c("CON"="#f6d55c", "MS"="#FFA080","ZVM"="#80A0FF")
metadata$TRT <- factor(metadata$TRT, levels = c("CON","MS","ZVM"))

## Unweighted UniFrac distance
Unwe_uni <- Unwe_uni_frac$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata)
Unwe_uni
min(Unwe_uni$PC1)
max(Unwe_uni$PC1)
min(Unwe_uni$PC2)
max(Unwe_uni$PC2)

a = ggplot(data = Unwe_uni,aes(x=PC1, y=PC2, color=TRT, fill=TRT)) +
  geom_point(shape = 21,alpha=1.0, size = 6, color = "black", stroke = 1.2) + #alpha controls transparency and helps when points are overlapping
  scale_fill_manual(values =  color) +
  scale_color_manual(values = color, guide = FALSE) +
  coord_cartesian(xlim = c(-0.35,0.35), ylim = c(-0.35,0.35)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2)+
  #theme_q2r()+
  stat_ellipse(geom="polygon", aes(fill=TRT), alpha=0.15, size=0.8, level = 0.79,type='norm')
a
a1 <- a+theme(axis.text.x = element_text(angle = 0, size = 26, family = "Times New Roman", colour = "black", face = "bold"), 
              axis.title.x = element_text(angle = 0, size = 28, family = "Times New Roman", colour = "black", face = "bold"),
              axis.title.y = element_text(size = 28, family = "Times New Roman", face = "bold"), legend.title = element_blank(),
              legend.text = element_text(size = 26, family = "Times New Roman", face = "bold", colour = "black"),
              axis.text.y = element_text(colour = "black", family = "Times New Roman", size = 26, face = "bold")) +
  labs(x = "PC1 (14.83%)", y = "PC2 (13.82%)", fill = "", face = 2) +
  theme(legend.position = "bottom")
a1
ggsave("230912_Soil_Unweighted_UniFrac.png", height=8, width=8, dpi = 300)

## Weighted UniFrac distance
we_uni <- we_uni_frac$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata)
we_uni
min(we_uni$PC1)
max(we_uni$PC1)
min(we_uni$PC2)
max(we_uni$PC2)

a = ggplot(data = we_uni,aes(x=PC1, y=PC2, color=TRT, fill=TRT)) +
  geom_point(shape = 21,alpha=1.0, size = 6, color = "black", stroke = 1.2) + #alpha controls transparency and helps when points are overlapping
  scale_fill_manual(values =  color) +
  scale_color_manual(values = color, guide = FALSE) +
  coord_cartesian(xlim = c(-0.1,0.1), ylim = c(-0.1,0.1)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  #theme_q2r()+
  stat_ellipse(geom="polygon", aes(fill=TRT), alpha=0.15, size=0.8, level = 0.80,type='norm')
a
a1 <- a+theme(axis.text.x = element_text(angle = 0, size = 26, family = "Times New Roman", colour = "black", face = "bold"), 
              axis.title.x = element_text(angle = 0, size = 28, family = "Times New Roman", colour = "black", face = "bold"),
              axis.title.y = element_text(size = 28, family = "Times New Roman", face = "bold"), legend.title = element_blank(),
              legend.text = element_text(size = 26, family = "Times New Roman", face = "bold", colour = "black"),
              axis.text.y = element_text(colour = "black", family = "Times New Roman", size = 26, face = "bold")) +
  labs(x = "PC1 (47.76%)", y = "PC2 (21.63%)", fill = "", face = 2) +
  theme(legend.position = "bottom")
a1
ggsave("230912_Soil_Weighted_UniFrac.png", height=8, width=8, dpi = 300)


# Taxa analysis
# https://github.com/jbisanz/qiime2R
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")
install.packages("ggplot2")
install.packages("reshape2")
install.packages('extrafont')
library(qiime2R)
library(dplyr)
library(tidyverse)
library("data.table")
# https://jkzorz.github.io/2019/06/05/stacked-bar-plots.html
library(ggplot2) # for plotting
library(reshape2) # for data manipulation
library(stringr)
# Change font
library(extrafont)
font_import()
fonts()
loadfonts(device="win")

setwd("D:/220729_Metagenomics_R_script/220808_R/230911_Hyunho_Soil_Metal/230918_Taxanomy_analysis")
list.files()
# How to use qza as a dataframe
SVs <- read_qza("table_filtered.qza")$data
SVs <- as.data.frame(SVs)
SVs <- tibble::rownames_to_column(SVs, "ASV")
metadata<-read_q2metadata("230911_Meta_data.txt")
taxonomy <- read_qza("taxonomy.qza")$data %>% parse_taxonomy(trim_extra = FALSE)
taxonomy <- as.data.frame(taxonomy)
View(taxonomy)
taxonomy[is.na(taxonomy)] <- "Unassigned"
taxonomy <- tibble::rownames_to_column(taxonomy, "ASV")
## Remove several ASVs removed from table_filtered.qza
taxa_SVs <- merge(SVs, taxonomy, by="ASV")
taxa_SVs
taxonomy <- taxa_SVs[c("ASV","Kingdom","Phylum","Class","Order","Family","Genus","Species")]
taxonomy <- as.data.frame(taxonomy)
## Change ASV to rownames from SVs and taxonomy
SVs_1 <- SVs[,-1]
rownames(SVs_1) <- SVs[,1]
SVs_1
taxonomy_1 <- taxonomy[,-1]
rownames(taxonomy_1) <- taxonomy[,1]
taxonomy_1
nrow(SVs_1)
nrow(taxonomy_1)

?str_split
# Phylum list
## Split the domain name by ";"
# taxasums_phylum$ASV <- str_split(taxasums_phylum$ASV,"; ", simplify = T)
taxasums_phylum <- summarize_taxa(SVs_1, taxonomy_1)$Phylum
taxasums_phylum <- tibble::rownames_to_column(taxasums_phylum, "ASV")
taxasums_phylum_data <- taxasums_phylum[,-1]
taxasums_phylum_name <- as.data.frame(taxasums_phylum[,1])
colnames(taxasums_phylum_name) <- c("Phylum")
taxasums_phylum <- cbind(taxasums_phylum_name,taxasums_phylum_data)
taxasums_phylum
nrow(taxasums_phylum)
write.table(taxasums_phylum,"230918_L2_Abs_all.txt", sep = "\t", row.names = FALSE)
## Occurrence based cut-off (30%)
n=nrow(metadata)
occurrence <- n-n*0.30
occurrence
taxasums_phylum_save <- taxasums_phylum[rowSums(taxasums_phylum == 0) <= occurrence, ]
nrow(taxasums_phylum_save)
## Raw-data for Ancon-BC analysis
write.table(taxasums_phylum_save, "230918_L2_Abs_occu_30.txt", sep = "\t", row.names = FALSE)

## Make relative abundance data: Taxa plot
taxasums_phylum <- as.data.frame(taxasums_phylum)
taxasums_phylum_1 <- taxasums_phylum[,-1]
rownames(taxasums_phylum_1) <- taxasums_phylum[,1]
taxasums_phylum_1
a <- t(taxasums_phylum_1)
View(a)
Phylum <- prop.table(a, margin = 1)
Phylum <- t(Phylum)
Phylum <- as.data.frame(Phylum)
Phylum_1 <- Phylum*100
sum(Phylum_1$CON_1)
View(Phylum_1)
nrow(Phylum_1)
## Occurrence based cut-off (30%)
B <- Phylum_1[rowSums(Phylum_1 == 0) <= occurrence, ]
nrow(B)
B <- rownames_to_column(B)
B
colnames(B)[1] <- "Phylum"
## Relative abundance data: non-low occurrence sequence variant
write.table(B,"230918_L2_RA_ocu_30.txt", sep = "\t", row.names = FALSE)

# Remove low RA taxa for bar plot
## Average based cut-off for barplot (0.5%)
colnames(B)[1] <- "rowname"
B <- column_to_rownames(B)
B <- B[(rowMeans(B) >= 1), ]
## Set the taxa by descending one
B <- B[order(rowMeans(B), decreasing = T), ]
B
nrow(B)

# If you need to change taxa name, then use this commend.
## In excel, just change taxa name.
B <- rownames_to_column(B)
B
colnames(B)[1] <- "Phylum"
write.table(B, "230918_L2_RA_occu30_RA_1.txt", sep = "\t", row.names = FALSE)

## After changing the taxa name, read the modified table using below code.
B <- read.table("230918_L2_RA_occu30_RA_1_Re.txt", header = TRUE, sep = "\t")
B
colnames(B)[1] <- "rowname"
B <- column_to_rownames(B)
Others <- as.data.frame(100-colSums(B))
colnames(Others) <- "Others"
Others <- t(Others)
B <- rbind(B,Others)
B

# Make input file for bar plot
B <- t(B)
B <- as.data.frame(B)
B <- tibble::rownames_to_column(B, "SampleID")
Treatment <- metadata[,1:2]
Phylum_2 <- merge(Treatment, B, by= "SampleID")
View(Phylum_2)

# Taxa bar plot at phylum level
# Change the format of table for plotting
AM = melt(Phylum_2, id = c("TRT","SampleID"))
AM
AM$Treatment = factor(AM$TRT, levels=c("CON","MS","ZVM"))
unique(AM$variable)
AM$Sample_ID = as.character(AM$SampleID)

colours_Dark2 = c("#1b9e77","#7570b3","#d95f02","#ff686b","#66a61e","#ffee93","#a6761d","#79addc","#b8e0d2","#ffee93","#ff686b",
                  "#1b9e77","#7570b3","#d95f02","#ff686b","#66a61e","#ffee93","#a6761d","#79addc","#b8e0d2","#ffee93","#ff686b",
                  "#1b9e77","#7570b3","#d95f02","#ff686b","#66a61e","#ffee93","#a6761d","#79addc","#b8e0d2","#ffee93","#ff686b",
                  "#1b9e77","#7570b3","#d95f02","#ff686b","#66a61e","#ffee93","#a6761d","#79addc","#b8e0d2","#ffee93","#ff686b",
                  "#1b9e77","#7570b3","#d95f02","#ff686b","#66a61e","#ffee93","#a6761d","#79addc","#b8e0d2","#ffee93","#ff686b",
                  "#1b9e77","#7570b3","#d95f02","#ff686b","#66a61e","#ffee93","#a6761d","#79addc","#b8e0d2","#ffee93","#ff686b")
Color <- unique(AM$variable)
Color_1 <- ifelse(Color == "Others", "grey40", colours_Dark2)
Color_1

# Make the plot
## scale_fill_brewer(palette = "Accent", "Pastel1", "Set1")
### http://rstudio-pubs-static.s3.amazonaws.com/5312_98fc1aba2d5740dd849a5ab797cc2c8d.html
mx <-  ggplot(AM, aes(x = SampleID, fill = variable, y = value))+
  geom_bar(stat = "identity", colour = "black", position="stack")+
  theme(axis.text.x = element_text(angle = 90, size = 22, family = "Times New Roman", 
                                   colour = "black", vjust = 0.5, hjust = 1.0, face = "bold"),
        axis.title.y = element_text(size = 24, family = "Times New Roman", face = "bold"),
        axis.title.x = element_text(size = 22, family = "Times New Roman", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 22, family = "Times New Roman", face = "bold", colour = "black"),
        axis.text.y = element_text(colour = "black", family = "Times New Roman", size = 22, face = "bold"),
        strip.text.x = element_text(size = 24, family = "Times New Roman", face = "bold"),
        plot.title = (element_text(size = 26, family = "Times New Roman", face = "bold", hjust = -0.08)),
        legend.position = "bottom")+
  labs(title = "(A) Phylum level", x = "Sample ID", y = "Relative Abundance (%)", fill = "")+
  scale_fill_manual(values = Color_1) +
  #scale_fill_manual(values = ifelse(AM$variable == "Others", "grey80", colours_Dark2)) +
  facet_grid(.~Treatment, scales = "free")+
  guides(fill=guide_legend(ncol = 4))
mx
ggsave('230918_L2_Taxa_bar_occu_30_ra_1.png', dpi = 600, height = 10, width = 15, unit = 'in')


df=ddply(AM, c("TRT","variable"), summarise, mean=mean(value), sd=sd(value), n= length(value), se=sd/sqrt(n))
df
colnames(df)[3]="Value"
df$Sample <- as.factor(df$Sample)
df$TRT <- factor(df$TRT, levels=c("CON", "MS", "ZVM"))
df

mx <-  ggplot(df, aes(x = TRT, fill = variable, y = Value))+
  geom_bar(stat = "identity", colour = "black", position="stack")+
  theme(axis.text.x = element_text(angle = 0, size = 22, family = "Times New Roman", 
                                   colour = "black", vjust = 0.5, hjust = 0.5, face = "bold"),
        axis.title.y = element_text(size = 24, family = "Times New Roman", face = "bold"),
        axis.title.x = element_text(size = 22, family = "Times New Roman", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 22, family = "Times New Roman", face = "bold", colour = "black"),
        axis.text.y = element_text(colour = "black", family = "Times New Roman", size = 22, face = "bold"),
        strip.text.x = element_text(size = 24, family = "Times New Roman", face = "bold"),
        plot.title = (element_text(size = 26, family = "Times New Roman", face = "bold", hjust = -0.75)),
        legend.position = "right")+
  labs(title = "(A) Phylum level", x = "", y = "Relative Abundance (%)", fill = "")+
  scale_fill_manual(values = Color_1) +
  #scale_fill_manual(values = ifelse(AM$variable == "Others", "grey80", colours_Dark2)) +
  #facet_grid(.~Treatment, scales = "free")+
  guides(fill=guide_legend(ncol = 1))
mx
ggsave('230918_L2_Taxa_bar_occu_30_ra_1_Grouped.png', dpi = 600, height = 8, width = 8, unit = 'in')


# Genus list
taxasums_Genus <- summarize_taxa(SVs_1, taxonomy_1)$Genus
taxasums_Genus <- tibble::rownames_to_column(taxasums_Genus, "ASV")
taxasums_Genus_data <- taxasums_Genus[,-1]
taxasums_Genus_name <- as.data.frame(taxasums_Genus[,1])
colnames(taxasums_Genus_name) <- c("Genus")
taxasums_Genus <- cbind(taxasums_Genus_name,taxasums_Genus_data)
taxasums_Genus
nrow(taxasums_Genus)
write.table(taxasums_Genus,"230918_L6_Abs_all.txt", sep = "\t", row.names = FALSE)
## Occurrence based cut-off (30%)
n=nrow(metadata)
occurrence <- n-n*0.30
occurrence
taxasums_Genus_save <- taxasums_Genus[rowSums(taxasums_Genus == 0) <= occurrence, ]
nrow(taxasums_Genus_save)
## Raw-data for Ancon-BC analysis
write.table(taxasums_Genus_save, "230918_L6_Abs_occu_30.txt", sep = "\t", row.names = FALSE)

## Make relative abundance data: Taxa plot
taxasums_Genus <- as.data.frame(taxasums_Genus)
taxasums_Genus_1 <- taxasums_Genus[,-1]
rownames(taxasums_Genus_1) <- taxasums_Genus[,1]
taxasums_Genus_1
a <- t(taxasums_Genus_1)
View(a)
Genus <- prop.table(a, margin = 1)
Genus <- t(Genus)
Genus <- as.data.frame(Genus)
Genus_1 <- Genus*100
sum(Genus_1$CON_1)
View(Genus_1)
nrow(Genus_1)
## Occurrence based cut-off (30%)
B <- Genus_1[rowSums(Genus_1 == 0) <= occurrence, ]
nrow(B)
B <- rownames_to_column(B)
B
colnames(B)[1] <- "Genus"
## Relative abundance data: non-low occurrence sequence variant
write.table(B,"230918_L6_RA_ocu_30.txt", sep = "\t", row.names = FALSE)

# Remove low RA taxa for bar plot
## Average based cut-off for barplot (0.5%)
colnames(B)[1] <- "rowname"
B <- column_to_rownames(B)
B <- B[(rowMeans(B) >= 1), ]
## Set the taxa by descending one
B <- B[order(rowMeans(B), decreasing = T), ]
B
nrow(B)

# If you need to change taxa name, then use this commend.
## In excel, just change taxa name.
B <- rownames_to_column(B)
B
colnames(B)[1] <- "Genus"
write.table(B, "230918_L6_RA_occu30_RA_1.txt", sep = "\t", row.names = FALSE)

## After changing the taxa name, read the modified table using below code.
B <- read.table("230918_L6_RA_occu30_RA_1_Re.txt", header = TRUE, sep = "\t")
B
colnames(B)[1] <- "rowname"
B <- column_to_rownames(B)
Others <- as.data.frame(100-colSums(B))
colnames(Others) <- "Others"
Others <- t(Others)
B <- rbind(B,Others)
B

# Make input file for bar plot
B <- t(B)
B <- as.data.frame(B)
B <- tibble::rownames_to_column(B, "SampleID")
Treatment <- metadata[,1:2]
Genus_2 <- merge(Treatment, B, by= "SampleID")
View(Genus_2)

# Taxa bar plot at Genus level
# Change the format of table for plotting
AM = melt(Genus_2, id = c("TRT","SampleID"))
AM
AM$Treatment = factor(AM$TRT, levels=c("CON","MS","ZVM"))
unique(AM$variable)
AM$Sample_ID = as.character(AM$SampleID)

colours_Dark2 = c("#1b9e77","#7570b3","#d95f02","#ff686b","#66a61e","#ffee93","#a6761d","#79addc","#b8e0d2","#ffee93","#ff686b",
                  "#1b9e77","#7570b3","#d95f02","#ff686b","#66a61e","#ffee93","#a6761d","#79addc","#b8e0d2","#ffee93","#ff686b",
                  "#1b9e77","#7570b3","#d95f02","#ff686b","#66a61e","#ffee93","#a6761d","#79addc","#b8e0d2","#ffee93","#ff686b",
                  "#1b9e77","#7570b3","#d95f02","#ff686b","#66a61e","#ffee93","#a6761d","#79addc","#b8e0d2","#ffee93","#ff686b",
                  "#1b9e77","#7570b3","#d95f02","#ff686b","#66a61e","#ffee93","#a6761d","#79addc","#b8e0d2","#ffee93","#ff686b",
                  "#1b9e77","#7570b3","#d95f02","#ff686b","#66a61e","#ffee93","#a6761d","#79addc","#b8e0d2","#ffee93","#ff686b")
Color <- unique(AM$variable)
Color_1 <- ifelse(Color == "Others", "grey40", colours_Dark2)
Color_1

# Make the plot
## scale_fill_brewer(palette = "Accent", "Pastel1", "Set1")
### http://rstudio-pubs-static.s3.amazonaws.com/5312_98fc1aba2d5740dd849a5ab797cc2c8d.html
mx <-  ggplot(AM, aes(x = SampleID, fill = variable, y = value))+
  geom_bar(stat = "identity", colour = "black", position="stack")+
  theme(axis.text.x = element_text(angle = 90, size = 22, family = "Times New Roman", 
                                   colour = "black", vjust = 0.5, hjust = 1.0, face = "bold"),
        axis.title.y = element_text(size = 24, family = "Times New Roman", face = "bold"),
        axis.title.x = element_text(size = 22, family = "Times New Roman", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 22, family = "Times New Roman", face = "bold", colour = "black"),
        axis.text.y = element_text(colour = "black", family = "Times New Roman", size = 22, face = "bold"),
        strip.text.x = element_text(size = 24, family = "Times New Roman", face = "bold"),
        plot.title = (element_text(size = 26, family = "Times New Roman", face = "bold", hjust = -0.08)),
        legend.position = "bottom")+
  labs(title = "(A) Phylum level", x = "Sample ID", y = "Relative Abundance (%)", fill = "")+
  scale_fill_manual(values = Color_1) +
  #scale_fill_manual(values = ifelse(AM$variable == "Others", "grey80", colours_Dark2)) +
  facet_grid(.~Treatment, scales = "free")+
  guides(fill=guide_legend(ncol = 4))
mx
ggsave('230918_L6_Taxa_bar_occu_30_ra_0.5.png', dpi = 600, height = 10, width = 15, unit = 'in')

df=ddply(AM, c("TRT","variable"), summarise, mean=mean(value), sd=sd(value), n= length(value), se=sd/sqrt(n))
df
colnames(df)[3]="Value"
df$TRT <- factor(df$TRT, levels=c("CON", "MS", "ZVM"))
df
label <- unique(df$variable)

mx <-  ggplot(df, aes(x = TRT, fill = variable, y = Value))+
  geom_bar(stat = "identity", colour = "black", position="stack")+
  theme(axis.text.x = element_text(angle = 0, size = 22, family = "Times New Roman", 
                                   colour = "black", vjust = 0.5, hjust = 0.5, face = "bold"),
        axis.title.y = element_text(size = 24, family = "Times New Roman", face = "bold"),
        axis.title.x = element_text(size = 22, family = "Times New Roman", face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 22, family = "Times New Roman", face = "bold", colour = "black"),
        axis.text.y = element_text(colour = "black", family = "Times New Roman", size = 22, face = "bold"),
        strip.text.x = element_text(size = 24, family = "Times New Roman", face = "bold"),
        plot.title = (element_text(size = 26, family = "Times New Roman", face = "bold", hjust = -0.63)),
        legend.position = "right")+
  labs(title = "(B) Genus level", x = "", y = "Relative Abundance (%)", fill = "")+
  scale_fill_manual(values = Color_1,
                    breaks = label,
                    labels = c(expression(bold("UCG_Gemmatimonadaceae")),
                               expression(bolditalic("Sphingomonas")),
                               expression(bolditalic("Bacillus")),
                               expression(bold("UCG_Chitinophagaceae")),
                               expression(bold("UG_Micrococcaceae")),
                               expression(bolditalic("Clostridium sensu stricto 1")),
                               expression(bold("JG30-KF-CM45")),
                               expression(bold("KD4-96")),
                               expression(bold("UCG_Microscillaceae")),
                               expression(bolditalic("Luteimonas")),
                               expression(bold("UCG_Vicinamibacterales")),
                               expression(bold("UG_Microbacteriaceae")),
                               expression(bolditalic("Longimicrobiaceae")),
                               expression(bolditalic("Vicinamibacteraceae")),
                               expression(bold("UG_Intrasporangiaceae")),
                               expression(bold("SC-I-84")),
                               expression(bolditalic("Lysobacter")),
                               expression(bolditalic("Bryobacter")),
                               expression(bold("Others"))))+
  guides(fill=guide_legend(ncol = 1))+
  theme(legend.text.align = 0)
mx
ggsave('230918_L6_Taxa_bar_occu_30_ra_1_Grouped.png', dpi = 600, height = 8, width = 9.1, unit = 'in')


# Venn diagram
# Install the packages
install.packages("VennDiagram")
library(VennDiagram)
library(ggplot2)
?venn.diagram

## Phylum level
### Before making the venn diagram, make input data (average by group).
### To make input file use this file: *_L2_Abs_all.txt (already produced)
### Remove Unassigned!

### Input data
# https://rdrr.io/github/yanlinlin82/ggvenn/man/geom_venn.html
## Raw data
setwd("D:/220729_Metagenomics_R_script/220808_R/230911_Hyunho_Soil_Metal/230918_Venn_diagram/")
A <- read.table("230918_Venn_L2.txt", sep = "\t", header = TRUE)
A
A1 <- A[,-1]
rownames(A1) <- A[,1]
colnames(A1) <- c("CON","MS","ZVM")
A1

## Treatment and Control should be present in Column
CON <- as.factor(rownames(A1)[A1[,"CON"] > 0])
MS <- as.factor(rownames(A1)[A1[,"MS"] > 0])
ZVM <- as.factor(rownames(A1)[A1[,"ZVM"] > 0])

# use list as input
A <-list('CON'=CON,'MS'=MS, 'ZVM'=ZVM)

## Check the image before save
venn.diagram(A, category.names = c("CON", "MS", "ZVM"),
             fill = c("#f6d55c", "#FFA080", "#80A0FF"), # color
             alpha = 0.8,
             col = "grey20", # edge color
             cex = 3, # size of the number
             lty = 2, # line type: 1. solid; 2. dashed; 3. dot; 4. duble dash; 5. long dash
             lwd = 3, # line tickness
             fontface = "bold", # face: plain, bold, italic
             fontfamily ="Times New Roman",
             cat.cex = 4, # size of the treatment label
             cat.col = c("black","black","black"),
             cat.alpha = 0.1,
             cat.fontface = "bold", # fontface of the label
             cat.fontfamily="Times New Roman",
             #cat.pos = c(-25,25,0), # label position (rotation)
             cat.dist = c(0.08, 0.08, -0.45), # label position (vertical)
             rotation.degree = 0, # set the rotation degree
             margin = 0.05, # set the figure size
             height = 30, width = 30, resolution = 600, imagetype = "png",
             units = "cm", disable.logging = TRUE,
             filename = "230918_Venn_L2_HKIM.png")

# Genus
A <- read.table("230918_Venn_L6.txt", sep = "\t", header = TRUE)
A
A1 <- A[,-1]
rownames(A1) <- A[,1]
colnames(A1) <- c("CON","MS","ZVM")
A1

## Treatment and Control should be present in Column
CON <- as.factor(rownames(A1)[A1[,"CON"] > 0])
MS <- as.factor(rownames(A1)[A1[,"MS"] > 0])
ZVM <- as.factor(rownames(A1)[A1[,"ZVM"] > 0])

# use list as input
A <-list('CON'=CON,'MS'=MS, 'ZVM'=ZVM)

## Check the image before save
venn.diagram(A, category.names = c("CON", "MS", "ZVM"),
             fill = c("#f6d55c", "#FFA080", "#80A0FF"), # color
             alpha = 0.8,
             col = "grey20", # edge color
             cex = 3, # size of the number
             lty = 2, # line type: 1. solid; 2. dashed; 3. dot; 4. duble dash; 5. long dash
             lwd = 3, # line tickness
             fontface = "bold", # face: plain, bold, italic
             fontfamily ="Times New Roman",
             cat.cex = 4, # size of the treatment label
             cat.col = c("black","black","black"),
             cat.alpha = 0.1,
             cat.fontface = "bold", # fontface of the label
             cat.fontfamily="Times New Roman",
             #cat.pos = c(-25,25,0), # label position (rotation)
             cat.dist = c(0.06, 0.06, 0.04), # label position (vertical)
             rotation.degree = 0, # set the rotation degree
             margin = 0.05, # set the figure size
             height = 30, width = 30, resolution = 600, imagetype = "png",
             units = "cm", disable.logging = TRUE,
             filename = "230918_Venn_L6_HKIM.png")

# Stat
## ANCOMBC2
library(ANCOMBC)
?ancombc2
setwd("D:/220729_Metagenomics_R_script/220808_R/230911_Hyunho_Soil_Metal/ANCOMBC2/230919_L2/")
### To analyze ancom-bc, we need phyloseq file.
### To make phyloseq file, we need three input files: 
### (1) Absolute abundance file, (2) Taxonomy file, and (3) Metadata
library("phyloseq")
# L2
otu_mat <- read.table("230919_Input_L2.txt", sep = "\t", header = TRUE)
tax_mat <- read.table("230919_L2_Taxonomy.txt", sep = "\t", header = TRUE)
samples_df <- read.table("230919_Metadata.txt", sep = "\t", header = TRUE)

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("ASV")
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample") 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

### Make phyloseq file
carbom <- phyloseq(OTU, TAX, samples)
sample_names(carbom)
rank_names(carbom)
sample_variables(carbom)

?ancombc
# Reference: HF-Low
set.seed(1234)
out = ancombc2(data = carbom, assay_name = NULL,
               tax_level = "Phylum", fix_formula = "TRT", 
               rand_formula = NULL, p_adj_method = "holm", 
               group = "TRT", prv_cut = 0.30,
               global = TRUE, pairwise = TRUE, 
               iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE), 
               em_control = list(tol = 1e-5, max_iter = 100),
               lme_control = NULL,
               mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
               struc_zero = FALSE, 
               alpha = 0.05, n_cl = 1, verbose = TRUE)
res_pair <- out$res_pair
res_global <- out$res_global
write.table(res_pair, "230919_L2_Ancombc_Pairwise_Soil.txt", sep = "\t")
write.table(res_global, "230919_L2_Ancombc_Global_Soil.txt", sep = "\t")

# L6
setwd("D:/220729_Metagenomics_R_script/220808_R/230911_Hyunho_Soil_Metal/ANCOMBC2/230919_L6/")
otu_mat <- read.table("230919_Input_L6.txt", sep = "\t", header = TRUE)
tax_mat <- read.table("230919_L6_Taxonomy.txt", sep = "\t", header = TRUE)
samples_df <- read.table("230919_Metadata.txt", sep = "\t", header = TRUE)

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("ASV")
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample") 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

### Make phyloseq file
carbom <- phyloseq(OTU, TAX, samples)
sample_names(carbom)
rank_names(carbom)
sample_variables(carbom)

?ancombc
set.seed(1234)
out = ancombc2(data = carbom, assay_name = NULL,
               tax_level = "Genus", fix_formula = "TRT", 
               rand_formula = NULL, p_adj_method = "holm", 
               group = "TRT", prv_cut = 0.30,
               global = TRUE, pairwise = TRUE, 
               iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE), 
               em_control = list(tol = 1e-5, max_iter = 100),
               lme_control = NULL,
               mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
               struc_zero = FALSE, 
               alpha = 0.05, n_cl = 1, verbose = TRUE)
res_pair <- out$res_pair
res_global <- out$res_global
write.table(res_pair, "230919_L6_Ancombc_Pairwise_Soil.txt", sep = "\t")
write.table(res_global, "230919_L6_Ancombc_Global_Soil.txt", sep = "\t")


# Stat_CON vs. MS
## ANCOMBC2
library(ANCOMBC)
?ancombc2
setwd("D:/220729_Metagenomics_R_script/220808_R/230911_Hyunho_Soil_Metal/ANCOMBC2/230919_L2/230926_CON_MS/")
### To analyze ancom-bc, we need phyloseq file.
### To make phyloseq file, we need three input files: 
### (1) Absolute abundance file, (2) Taxonomy file, and (3) Metadata
library("phyloseq")
# L2
otu_mat <- read.table("230919_Input_L2.txt", sep = "\t", header = TRUE)
tax_mat <- read.table("230919_L2_Taxonomy.txt", sep = "\t", header = TRUE)
samples_df <- read.table("230919_Metadata.txt", sep = "\t", header = TRUE)

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("ASV")
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample") 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

### Make phyloseq file
carbom <- phyloseq(OTU, TAX, samples)
sample_names(carbom)
rank_names(carbom)
sample_variables(carbom)

??ancombc
# Reference: HF-Low
set.seed(1234)
out = ancombc2(data = carbom, assay_name = NULL,
               tax_level = "Phylum", fix_formula = "TRT", 
               rand_formula = NULL, p_adj_method = "holm", 
               group = "TRT", prv_cut = 0.30,
               global = FALSE, pairwise = FALSE, 
               iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE), 
               em_control = list(tol = 1e-5, max_iter = 100),
               lme_control = NULL,
               mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
               struc_zero = FALSE, 
               alpha = 0.05, n_cl = 1, verbose = TRUE)
res <- out$res
write.table(res, "230919_L2_Ancombc_Pairwise_Soil_CON_MS.txt", sep = "\t")

# L6
setwd("D:/220729_Metagenomics_R_script/220808_R/230911_Hyunho_Soil_Metal/ANCOMBC2/230919_L6/230926_CON_VS_MS/")
otu_mat <- read.table("230919_Input_L6.txt", sep = "\t", header = TRUE)
tax_mat <- read.table("230919_L6_Taxonomy.txt", sep = "\t", header = TRUE)
samples_df <- read.table("230919_Metadata.txt", sep = "\t", header = TRUE)

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("ASV")
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample") 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

### Make phyloseq file
carbom <- phyloseq(OTU, TAX, samples)
sample_names(carbom)
rank_names(carbom)
sample_variables(carbom)

?ancombc
set.seed(1234)
out = ancombc2(data = carbom, assay_name = NULL,
               tax_level = "Genus", fix_formula = "TRT", 
               rand_formula = NULL, p_adj_method = "holm", 
               group = "TRT", prv_cut = 0.30,
               global = FALSE, pairwise = FALSE, 
               iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE), 
               em_control = list(tol = 1e-5, max_iter = 100),
               lme_control = NULL,
               mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
               struc_zero = FALSE, 
               alpha = 0.05, n_cl = 1, verbose = TRUE)
res <- out$res
write.table(res, "230919_L6_Ancombc_Pairwise_Soil_CON_MS.txt", sep = "\t")

# KEGG_Modules
setwd("D:/220729_Metagenomics_R_script/220808_R/230911_Hyunho_Soil_Metal/ANCOMBC2/231010_KEGG_BRITE/")
### To analyze ancom-bc, we need phyloseq file.
### To make phyloseq file, we need three input files: 
### (1) Absolute abundance file, (2) Taxonomy file, and (3) Metadata
library("phyloseq")
library(reshape2)
library(dplyr)

# KEGG_Modules
otu_mat <- read.table("230919_Input_L2.txt", sep = "\t", header = TRUE)
tax_mat <- read.table("230919_L2_Taxonomy.txt", sep = "\t", header = TRUE)
samples_df <- read.table("230919_Metadata.txt", sep = "\t", header = TRUE)

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("ASV")
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample") 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

### Make phyloseq file
carbom <- phyloseq(OTU, TAX, samples)
sample_names(carbom)
rank_names(carbom)
sample_variables(carbom)

?ancombc
# Reference: HF-Low
set.seed(1234)
out = ancombc2(data = carbom, assay_name = NULL,
               tax_level = "Phylum", fix_formula = "TRT", 
               rand_formula = NULL, p_adj_method = "holm", 
               group = "TRT", prv_cut = 0,
               global = FALSE, pairwise = FALSE, 
               iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE), 
               em_control = list(tol = 1e-5, max_iter = 100),
               lme_control = NULL,
               mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
               struc_zero = FALSE, 
               alpha = 0.05, n_cl = 1, verbose = TRUE)
res <- out$res
write.table(res, "231010_KEGG_Modules_Ancombc_Pairwise_Soil_CON_MS.txt", sep = "\t")

# KEGG_Pathway
setwd("D:/220729_Metagenomics_R_script/220808_R/230911_Hyunho_Soil_Metal/ANCOMBC2/231010_KEGG_PATHWAY/")
### To analyze ancom-bc, we need phyloseq file.
### To make phyloseq file, we need three input files: 
### (1) Absolute abundance file, (2) Taxonomy file, and (3) Metadata
library("phyloseq")
library(reshape2)
library(dplyr)

# KEGG_Pathway
otu_mat <- read.table("230919_Input_L2.txt", sep = "\t", header = TRUE)
tax_mat <- read.table("230919_L2_Taxonomy.txt", sep = "\t", header = TRUE)
samples_df <- read.table("230919_Metadata.txt", sep = "\t", header = TRUE)

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("ASV")
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample") 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

### Make phyloseq file
carbom <- phyloseq(OTU, TAX, samples)
sample_names(carbom)
rank_names(carbom)
sample_variables(carbom)

?ancombc
# Reference: HF-Low
set.seed(1234)
out = ancombc2(data = carbom, assay_name = NULL,
               tax_level = "Phylum", fix_formula = "TRT", 
               rand_formula = NULL, p_adj_method = "holm", 
               group = "TRT", prv_cut = 0,
               global = FALSE, pairwise = FALSE, 
               iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE), 
               em_control = list(tol = 1e-5, max_iter = 100),
               lme_control = NULL,
               mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
               struc_zero = FALSE, 
               alpha = 0.05, n_cl = 1, verbose = TRUE)
res <- out$res
write.table(res, "231010_KEGG_Pathway_Ancombc_Pairwise_Soil_CON_MS.txt", sep = "\t")

# KEGG_Orthology_PCA
# PCA analysis with PERMANOVA
## https://rpubs.com/collnell/manova
## https://github.com/pmartinezarbizu/pairwiseAdonis
install.packages('vegan')
install.packages("ggplot2")
install.packages("ggfortify")
devtools::install_github("cmartin/ggConvexHull")
install.packages('devtools')
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(devtools)
library("pairwiseAdonis")
library(ggplot2)
library(ggfortify)
library(ggConvexHull)
library(vegan) ##Community ecology: ordination, disversity & dissimilarities
# PICRUSt2 analysis
# PICRUSt2_KO_metagenome_out/pred_metagenome_contrib.tsv
library(dplyr)
library(reshape2)
library(qiime2R)
setwd("D:/220729_Metagenomics_R_script/220808_R/230911_Hyunho_Soil_Metal/230920_KEGG/230920_PCA/")
metadata<-read_q2metadata("230911_Meta_data.txt")
metadata
# Using pred_metagenome_unstrat.tsv file
KO <- read.table("pred_metagenome_unstrat.tsv", sep = "\t", header = TRUE)
colnames(KO)[1] <- "SampleID"
KO_1 <- t(KO)
KO_1 <- as.data.frame(KO_1)
KO_2 <- KO_1[-1,]
KO_2 <- rownames_to_column(KO_2)
Colname <- KO_1[1,]
Colname <- rownames_to_column(Colname)
colnames(KO_2) <- Colname
View(KO_2)
KOs <- merge(metadata[,1:2],KO_2, by="SampleID")
View(KOs)

# PCA analysis
# set working directory to save output file!
KOs.matrix <- as.matrix(KOs[,-c(1:2)])
View(KOs.matrix)
KOs.prop <- decostand(KOs.matrix, method="total")
KOs.dist<-vegdist(KOs.prop, method='bray')
set.seed(36)
View(KOs)
KOs.div <- adonis2(KOs.dist~TRT, data=KOs, permutations = 9999, method="bray")
KOs.div
set.seed(36)
?pairwise.adonis2
KOs.pairwise <- pairwise.adonis2(KOs.dist~TRT, data=KOs, perm = 9999, p.adjust.methods = "bonferroni")
KOs.pairwise

# remove treatment line from A
A_dt <- KOs[,-c(1:2)]
View(A_dt)
# make treatment data from A
A_dt_group <- KOs[,2]
A_dt_group
# Do PCA analysis
# https://stats.stackexchange.com/questions/53/pca-on-correlation-or-covariance
# There are two types of standard for PCA analysis, Covariance matrix (scale = FALSE) and correlation matrix (scale = TRUE)
# For the normalized data, we have to use result from Covariance matrix
PCA_result <- prcomp(A_dt, center = T, scale. = FALSE)
PCA_result
# Print proportioin of variance / Check elbow point
plot(PCA_result, type = "l")
# Check the total proportion from PC1 to PC nubmer at elbow point
# General standard for elbow point = 90%
summary(PCA_result)
?stat_ellipse
g=autoplot(PCA_result, data = KOs, colour = 'TRT',shape ='TRT', fill='TRT',frame = FALSE, frame.type = 'norm', 
           size = 8, )+
  stat_ellipse(geom = "polygon", aes(fill = TRT), level = 0.95, alpha = 0.15)+
  ggtitle("")+
  theme(axis.text.x = element_text(angle = 0, size = 24, family = "Times New Roman", colour = "black", face = "bold"), 
        axis.title.x = element_text(angle = 0, size = 26, family = "Times New Roman", colour = "black", face = "bold"),
        axis.title.y = element_text(size = 26, family = "Times New Roman", face = "bold"), 
        legend.title = element_blank(),
        legend.text = element_text(size = 26, family = "Times New Roman", face = "bold", colour = "black"),
        axis.text.y = element_text(colour = "black", family = "Times New Roman", size = 24, face = "bold"),
        plot.title = element_text(size = 28, family = "Times New Roman", face = "bold"))+
  scale_shape_manual(breaks=c("CON", "MS", "ZVM"),
                     values=c("CON" = 21,"MS" = 21, "ZVM" = 21))+
  scale_fill_manual(breaks=c("CON", "MS", "ZVM"),
                    values=c("CON" = "#ffee93","MS" = "#b8e0d2", "ZVM" = "#7570b3"))  +
  scale_colour_manual(breaks=c("CON", "MS", "ZVM"),
                      values=c("CON" = "black","MS" = "black", "ZVM" = "black"))  +
  geom_hline(yintercept=0, linetype="dashed", color = "grey60", size = 0.8) +
  geom_vline(xintercept=0, linetype="dashed", colour= "grey60", size = 0.8)
g
g$layers[[1]]$aes_params$size <- 8.5
g
g1=g + theme(legend.position = "bottom")
g1 <- g1 + theme(plot.title = element_text(hjust = -0.25))
g1
ggsave("230920_PCA_KOs.png", width = 9, height = 9, units = "in", dpi = 600)


# KEGG Brite Module classification
library(dplyr)
library(reshape)
library(qiime2R)
library(plyr)
library(reshape2)
setwd("D:/220729_Metagenomics_R_script/220808_R/230911_Hyunho_Soil_Metal/230920_KEGG/231010_KEGG_BRITE")
# Input data
A <- read.table("231010_Input.txt", sep = "\t", header = TRUE)
A

# Metadata
metadata

label <- A[,1:6]
df <- A[,-c(2:6)]
Long = melt(df, id = c("Module"))
Long_merged <- merge(Long, label, by = "Module")
View(Long_merged)

# Brite_Level_1
Long_merged_L1 <- Long_merged[,c("Brite_L1", "variable", "value")]
L1 <- Long_merged_L1 %>% aggregate(value~Brite_L1+variable, FUN = sum)
L1_RA <- ddply(L1, .(variable), mutate, L1_RA = value / sum(value)*100)
L1_RA <- L1_RA[,c(1,2,4)]
colnames(L1_RA) <- c("Brite_L1","SampleID","L1_RA")
View(L1_RA)

# wide_form: General format
Wide_L1_RA <- dcast(L1_RA, Brite_L1 ~ SampleID, value.var = "L1_RA")
View(Wide_L1_RA)
write.table(Wide_L1_RA, "231010_KEGG_Brite_L1_RA_Wide.txt", sep = "\t")

metadata_1 <- metadata[,c(1:2)]
L1_RA <- merge(L1_RA, metadata_1, by = "SampleID")
L1_Stat=ddply(L1_RA, c("Brite_L1","TRT"), summarise, mean=mean(L1_RA), sd=sd(L1_RA), n= length(L1_RA), se=sd/sqrt(n))
L1_Stat # for making several plots.

# Brite_Level_2
Long_merged_L2 <- Long_merged[,c("Brite_L2", "variable", "value")]
L2 <- Long_merged_L2 %>% aggregate(value~Brite_L2+variable, FUN = sum)
L2_RA <- ddply(L2, .(variable), mutate, L2_RA = value / sum(value)*100)
L2_RA <- L2_RA[,c(1,2,4)]
colnames(L2_RA) <- c("Brite_L2","SampleID","L2_RA")
View(L2_RA)

# wide_form: General format
Wide_L2_RA <- dcast(L2_RA, Brite_L2 ~ SampleID, value.var = "L2_RA")
View(Wide_L2_RA)
write.table(Wide_L2_RA, "231010_KEGG_Brite_L2_RA_Wide.txt", sep = "\t")

metadata_1 <- metadata[,c(1:2)]
L2_RA <- merge(L2_RA, metadata_1, by = "SampleID")
L2_Stat=ddply(L2_RA, c("Brite_L2","TRT"), summarise, mean=mean(L2_RA), sd=sd(L2_RA), n= length(L2_RA), se=sd/sqrt(n))
L2_Stat # for making several plots.

# Brite_Level_3
Long_merged_L3 <- Long_merged[,c("Brite_L3", "variable", "value")]
L3 <- Long_merged_L3 %>% aggregate(value~Brite_L3+variable, FUN = sum)
L3_RA <- ddply(L3, .(variable), mutate, L3_RA = value / sum(value)*100)
L3_RA <- L3_RA[,c(1,2,4)]
colnames(L3_RA) <- c("Brite_L3","SampleID","L3_RA")
View(L3_RA)

# wide_form: General format
Wide_L3_RA <- dcast(L3_RA, Brite_L3 ~ SampleID, value.var = "L3_RA")
View(Wide_L3_RA)
write.table(Wide_L3_RA, "231010_KEGG_Brite_L3_RA_Wide.txt", sep = "\t")

metadata_1 <- metadata[,c(1:2)]
L3_RA <- merge(L3_RA, metadata_1, by = "SampleID")
L3_Stat=ddply(L3_RA, c("Brite_L3","TRT"), summarise, mean=mean(L3_RA), sd=sd(L3_RA), n= length(L3_RA), se=sd/sqrt(n))
L3_Stat # for making several plots.

# Brite_Level_4
Long_merged_L4 <- Long_merged[,c("Brite_L4", "variable", "value")]
L4 <- Long_merged_L4 %>% aggregate(value~Brite_L4+variable, FUN = sum)
L4_RA <- ddply(L4, .(variable), mutate, L4_RA = value / sum(value)*100)
L4_RA <- L4_RA[,c(1,2,4)]
colnames(L4_RA) <- c("Brite_L4","SampleID","L4_RA")
View(L4_RA)

# wide_form: General format
Wide_L4_RA <- dcast(L4_RA, Brite_L4 ~ SampleID, value.var = "L4_RA")
View(Wide_L4_RA)
write.table(Wide_L4_RA, "231010_KEGG_Brite_L4_RA_Wide.txt", sep = "\t")

metadata_1 <- metadata[,c(1:2)]
L4_RA <- merge(L4_RA, metadata_1, by = "SampleID")
L4_Stat=ddply(L4_RA, c("Brite_L4","TRT"), summarise, mean=mean(L4_RA), sd=sd(L4_RA), n= length(L4_RA), se=sd/sqrt(n))
L4_Stat # for making several plots.

L1_Stat
L2_Stat
L3_Stat
L4_Stat
write.table(L1_Stat,"231010_L1_Stat.txt", sep = "\t", row.names = FALSE)
write.table(L2_Stat,"231010_L2_Stat.txt", sep = "\t", row.names = FALSE)
write.table(L3_Stat,"231010_L3_Stat.txt", sep = "\t", row.names = FALSE)
write.table(L4_Stat,"231010_L4_Stat.txt", sep = "\t", row.names = FALSE)


# Correlation analysis between species and kegg modules
## https://microbiome.github.io/tutorials/Heatmap.html
###Installation of the microbiome package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")
library(microbiome)

# Check the input file
setwd("D:/220729_Metagenomics_R_script/220808_R/230911_Hyunho_Soil_Metal/231010_Correaltion/")
Microbes <- read.table("231010_Input_Microbes.txt", sep = "\t", header = TRUE)
KEGG_Module <- read.table("231010_Input_Soil.txt", sep = "\t", header = TRUE)
x <- Microbes
y <- KEGG_Module

correlations <- associate(x,y, method="spearman", mode="matrix", p.adj.threshold = 1.0, n.signif = 1,
                          p.adj.method = "none")
correlations
correlation.table <- associate(x,y, method="spearman", mode="table", p.adj.threshold = 1.0, n.signif = 1,
                               p.adj.method = "none")
View(correlation.table)

#This is to change correlation table into txt.file
write.table(correlation.table, file="Result_microbiome_corr_table_not_adjusted.txt", append= TRUE)

# Below 3 methods is for making correlation heatmap.  
# method_1 of microbiome package_using_association heatmaps
p <- heat(correlation.table, "X1","X2",fill = "Correlation",star = "p.adj",p.adj.threshold = 0.05)
print(p)

install.packages("ggplot2")
library(ggplot2)
library(dplyr)
theme_set(theme_bw(20))
df <- correlation.table
df <- df %>% mutate(label = case_when(p.adj > 0.05 ~ "",
                                      p.adj > 0.01 ~ "*",
                                      p.adj > 0.001 ~ "**",
                                      !is.na(p.adj) ~ "***",
                                      TRUE ~ NA_character_))
View(df)
write.table(df, file = "231010_Significant_Corr_Speaman_L2_L6.txt", sep = "\t")
?scale_fill_gradientn
p <- ggplot(df, aes(X1, X2, group=X2)) +
  geom_tile(aes(fill = Correlation)) +
  #geom_text(aes(label = round(df$Correlation, 1)), size = 2)+
  scale_fill_gradientn("Correlation", 
                       breaks = seq(from = -1, to = 1,  by = 0.25), 
                       colours = c("#3988C6", "#F0F2E4", "#f65A46"), 
                       limits = c(-1, 1), aesthetics = "fill",
                       space = "Lab",
                       name = "Spearman\ncorrelation") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 30, face="bold", family = "Times New Roman", angle = 90, vjust = 0.4, hjust = 1.0), 
        axis.text.y = element_text(colour = "black", face="bold", family = "Times New Roman", size = 30), 
        legend.text = element_text(size = 30, face="bold", family = "Times New Roman", colour ="black"), 
        legend.title = element_text(size = 30, face="bold", family = "Times New Roman"),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey85"),
        plot.title = element_text(colour = "black", size = 16, face="bold", family = "Times New Roman")) +
  geom_text(aes(x = X1, y = X2, label = label, vjust = 0.75, hjust = 0.5), col = "black", size = 12) +
  xlab("") + ylab("") +
  scale_x_discrete(limits = c("p_Bacteroidota",
                              "p_Hydrogenedentes",
                              "UG_Micrococcaceae",
                              "g_Luteimonas",
                              "g_Longimicrobiaceae",
                              "g_Gitt.GS.136",
                              "UG_Sphingomonadaceae",
                              "g_Skermanella",
                              "g_Nocardioides",
                              "UG_Gemmatimonadaceae",
                              "g_JG30.KF.CM66",
                              "g_Nannocystis",
                              "g_AKYG1722",
                              "g_AKAU4049",
                              "UCG_Anaerolineaceae",
                              "g_Nitrospira",
                              "g_Ureibacillus",
                              "g_Herminiimonas",
                              "UCG_Elsterales",
                              "g_Magnetospirillaceae",
                              "g_Pricia",
                              "g_B1.7BS",
                              "g_Afipia",
                              "UCG_Cyclobacteriaceae",
                              "g_Anaeromyxobacter",
                              "UCG_Rhodospirillales"),
                   labels = c(expression(bold("p_Bacteroidota")),
                              expression(bold("p_Hydrogenedentes")),
                              expression(bold("UG_Micrococcaceae")),
                              expression(bold("g")*bolditalic("Luteimonas")),
                              expression(bold("g")*bolditalic("Longimicrobiaceae")),
                              expression(bold("g_Gitt-GS-136")),
                              expression(bold("UG_Sphingomonadaceae")),
                              expression(bold("g")*bolditalic("Skermanella")),
                              expression(bold("g")*bolditalic("Nocardioides")),
                              expression(bold("UG_Gemmatimonadaceae")),
                              expression(bold("g_JG30-KF-CM66")),
                              expression(bold("g")*bolditalic("Nannocystis")),
                              expression(bold("g_AKYG1722")),
                              expression(bold("g_AKAU4049")),
                              expression(bold("UCG_Anaerolineaceae")),
                              expression(bold("g")*bolditalic("Nitrospira")),
                              expression(bold("g")*bolditalic("Ureibacillus")),
                              expression(bold("g")*bolditalic("Herminiimonas")),
                              expression(bold("UCG_Elsterales")),
                              expression(bold("g")*bolditalic("Magnetospirillaceae")),
                              expression(bold("g")*bolditalic("Pricia")),
                              expression(bold("g_B1-7BS")),
                              expression(bold("g")*bolditalic("Afipia")),
                              expression(bold("UCG_Cyclobacteriaceae")),
                              expression(bolditalic("Anaeromyxobacter")),
                              expression(bold("UCG_Rhodospirillales"))))+
  scale_y_discrete(limits = c("hao", "amoB","amoA", "nosZ","norB","nirS",
                              "nirK","narG","dsrA","Soil.pH","SO4","NO3.","NH4.","N2O"),
                   labels = c(expression(bolditalic("hao")),
                              expression(bolditalic("amoB")),
                              expression(bolditalic("amoA")),
                              expression(bolditalic("nosZ")),
                              expression(bolditalic("norB")),
                              expression(bolditalic("nirS")),
                              expression(bolditalic("nirK")),
                              expression(bolditalic("narG")),
                              expression(bolditalic("dsrA")),
                              expression(bold("Soil pH")),
                              expression(bold(SO[4]^{" 2-"})),
                              expression(bold(NO[3]^{" -"})),
                              expression(bold(NH[4]^{" +"})),
                              expression(bold(N[2]*O))))+
  theme(legend.key.width = unit(2, "cm"),
        legend.key.height = unit(2.2, "cm")) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))
p
ggsave("231010_Sig_Genera_Soil_properties_Spearman_L2_L6.png", width = 20, height = 16, dpi = 600)


# ANCOM-BC Plot
# Genus
## Input
setwd("D:/220729_Metagenomics_R_script/220808_R/230911_Hyunho_Soil_Metal/ANCOMBC2/231010_Plot/231010_L6/")
A <- read.table("231010_Plot_L6_Ancombc_Holm_Bonf_Sig.txt", sep = "\t", header = TRUE)

# CON vs. MS
CON.MS <- A %>% subset(A$Variable == "CON vs. MS")
CON.MS$col <- gsub("blue","#3988C6",CON.MS$col)
CON.MS$col <- gsub("red1","#f65A46",CON.MS$col)
name=c("1 : plane","2 : bold","3 : italic","4 : bold and italic")
g = ggplot(CON.MS, aes(x = LFC, y = Taxa, fill = Variable, color = Variable))+
  geom_errorbarh(aes(xmin=Lower,xmax=Upper),color="black", height = 0.1) +
  geom_point(shape = 18, size = 6.0, color = CON.MS$col)+
  scale_x_continuous(breaks = seq(-5,6.5,1))+
  geom_vline(xintercept = 0, color="grey55", size = 1.2, linetype = "dashed")+
  ggtitle("")+
  theme_grey()
g
g1 <- g + theme(legend.title = element_blank(), 
                legend.text = element_blank(),
                legend.position = "none",
                axis.title.y = element_blank(),
                axis.text.x = element_text(face = "bold",colour = "black", size = 15, family = "Times New Roman"), 
                axis.text.y = element_text(face = "bold",colour = "black", size = 15, family = "Times New Roman"), 
                axis.title.x = element_text(face = "bold", size = 16, colour = "black", family = "Times New Roman",
                                            margin = margin(t = 0, r = 0, b = 25, l = 0)),
                strip.text.x = element_text(face = "bold", size = 15, colour = "black", family = "Times New Roman"),
                plot.title = element_text(face = "bold", size = 18, colour = "black", family = "Times New Roman"),
                axis.line = element_line(colour = "black", size = 0.8))+ 
  labs(x= "", y = "", fill = "")+
  coord_cartesian(xlim = c(-2.5,6.25), ylim = c(1,24), clip = "off")+ # change ylim number matching the number of taxa
  annotate("text", x = 3.5, y = -1.8, label = "Increased in MS",
           family = "Times New Roman", fontface = 2, size = 5.0)+
  annotate("text", x = 0.0, y = -0.7, label = "Log fold change",
           family = "Times New Roman", fontface = 2, size = 5.5)+
  annotate("segment", x = 0.25, y = -1.3, xend = 6.5, yend = -1.3, 
           size = 1.4, lineend = "round", linejoin = "round",
           arrow = arrow(length = unit(0.15, "inches")))+
  annotate("text", x = -1.5, y = -1.8, label = "Increased in CON",
           family = "Times New Roman", fontface = 2, size = 5.0)+
  annotate("segment", x = -0.2, y = -1.3, xend = -3.0, yend = -1.3, 
           size = 1.4, lineend = "round", linejoin = "round",
           arrow = arrow(length = unit(0.15, "inches")))+
  scale_y_discrete(limits = c("UCG_Rhodospirillales","Anaeromyxobacter","UCG_Cyclobacteriaceae","Afipia",
                              "B1-7BS","Pricia","Magnetospirillaceae","UCG_Elsterales","Herminiimonas",
                              "Ureibacillus","Nitrospira","UCG_Anaerolineaceae","AKAU4049","AKYG1722",
                              "Nannocystis","JG30-KF-CM66","UG_Gemmatimonadaceae","Nocardioides",
                              "Skermanella","UG_Sphingomonadaceae","Gitt-GS-136",
                              "Longimicrobiaceae","Luteimonas","UG_Micrococcaceae"),
                   labels = c(expression(bold("UCG_Rhodospirillales")),
                              expression(bolditalic("Anaeromyxobacter")),
                              expression(bold("UCG_Cyclobacteriaceae")),
                              expression(bolditalic("Afipia")),
                              expression(bold("B1-7BS")),
                              expression(bolditalic("Pricia")),
                              expression(bolditalic("Magnetospirillaceae")),
                              expression(bold("UCG_Elsterales")),
                              expression(bolditalic("Herminiimonas")),
                              expression(bolditalic("Ureibacillus")),
                              expression(bolditalic("Nitrospira")),
                              expression(bold("UCG_Anaerolineaceae")),
                              expression(bold("AKAU4049")),
                              expression(bold("AKYG1722")),
                              expression(bolditalic("Nannocystis")),
                              expression(bold("JG30-KF-CM66")),
                              expression(bold("UG_Gemmatimonadaceae")),
                              expression(bolditalic("Nocardioides")),
                              expression(bolditalic("Skermanella")),
                              expression(bold("UG_Sphingomonadaceae")),
                              expression(bold("Gitt-GS-136")),
                              expression(bolditalic("Longimicrobiaceae")),
                              expression(bolditalic("Luteimonas")),
                              expression(bold("UG_Micrococcaceae"))))
g1

p_right <- ggplot(data = CON.MS, aes(x=0, y=Taxa, label=P)) +
  theme_void()+
  geom_text(aes(x = 0, y = Taxa, label = P),size = 5.2,
            hjust = 0.5, family = "Times New Roman", fontface="bold") +
  ggtitle(expression(atop(paste(bold("Adjusted")~bolditalic("P")))))+
  theme(plot.title = element_text(hjust = 0.5, vjust = -6, size = 15, face = 2, family = "Times New Roman"),
        legend.text = element_blank(),
        legend.position = "none")
p_right

library("patchwork")
layout <- c(
  area(t = 0, l = 0, b = 10, r = 6.9), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  area(t = 0, l = 6.0, b = 10, r = 8) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
)
# final plot arrangement
P_CON.MS <- g1 + theme(plot.title = element_text(hjust = -2.1)) + p_right + plot_layout(design = layout)
P_CON.MS
ggsave("231010_ANCOM-BC_L6_CON_MS.png",
       width = 10, height = 10, units = "in", dpi = 600)

# L6 Relative abundance bar plot
library("plyr")
library(reshape2)
B <-  read.table("231010_L6_bar_plot_Relative_abundance.txt", sep="\t", header = T)
# Change the format of table for plotting
colnames(B) <- c("Sample","CON","MS")
B
AM = melt(B, id = c("Sample"))
AM
AM$variable <- factor(AM$variable, levels = c("CON","MS"))
df=ddply(AM, c("Sample","variable"), summarise, mean=mean(value), sd=sd(value), n= length(value), se=sd/sqrt(n))
df
colnames(df)[3]="Value"
df
mx <-  ggplot(df, aes(Sample, Value, fill=variable))+
  theme_gray()

# font --> face = 1 ---
par(mfrow=c(2,2))
name=c("1 : plane","2 : bold","3 : italic","4 : bold and italic")
TRT_color <- c("CON"="#f6d55c", "MS"="#FFA080")
# phylum to genus
p1=mx+geom_bar(stat="identity", position=position_dodge(-0.8), width = 0.8)+
  coord_flip(ylim = c(0,3.5))+
  scale_x_discrete(limits = c("UCG_Rhodospirillales","Anaeromyxobacter","UCG_Cyclobacteriaceae","Afipia",
                              "B1-7BS","Pricia","Magnetospirillaceae","UCG_Elsterales","Herminiimonas",
                              "Ureibacillus","Nitrospira","UCG_Anaerolineaceae","AKAU4049","AKYG1722",
                              "Nannocystis","JG30-KF-CM66","UG_Gemmatimonadaceae","Nocardioides",
                              "Skermanella","UG_Sphingomonadaceae","Gitt-GS-136",
                              "Longimicrobiaceae","Luteimonas","UG_Micrococcaceae"),
                   labels = c(expression(bold("UCG_Rhodospirillales")),
                              expression(bolditalic("Anaeromyxobacter")),
                              expression(bold("UCG_Cyclobacteriaceae")),
                              expression(bolditalic("Afipia")),
                              expression(bold("B1-7BS")),
                              expression(bolditalic("Pricia")),
                              expression(bolditalic("Magnetospirillaceae")),
                              expression(bold("UCG_Elsterales")),
                              expression(bolditalic("Herminiimonas")),
                              expression(bolditalic("Ureibacillus")),
                              expression(bolditalic("Nitrospira")),
                              expression(bold("UCG_Anaerolineaceae")),
                              expression(bold("AKAU4049")),
                              expression(bold("AKYG1722")),
                              expression(bolditalic("Nannocystis")),
                              expression(bold("JG30-KF-CM66")),
                              expression(bold("UG_Gemmatimonadaceae")),
                              expression(bolditalic("Nocardioides")),
                              expression(bolditalic("Skermanella")),
                              expression(bold("UG_Sphingomonadaceae")),
                              expression(bold("Gitt-GS-136")),
                              expression(bolditalic("Longimicrobiaceae")),
                              expression(bolditalic("Luteimonas")),
                              expression(bold("UG_Micrococcaceae"))))+
  geom_errorbar(aes(ymin=Value-se, ymax=Value+se),
                position=position_dodge(-0.8), width=0.1)+
  ggtitle("")+
  theme(axis.text.x = element_text(angle = 0, size = 15, family = "Times New Roman", colour = "black", face = "bold"), 
        axis.title.x = element_text(size = 15, family = "Times New Roman", face = "bold"),
        axis.title.y = element_text(size = 15, family = "Times New Roman", face = "bold"), 
        legend.title = element_blank(),
        legend.text = element_blank(),
        axis.text.y = element_text(colour = "black", family = "Times New Roman", size = 15, face=2),
        plot.title = element_text(size = 18, family = "Times New Roman", face = "bold"),
        text = element_text(size = 15, family = "Times New Roman", face = "bold", colour = "black",
                            vjust = 0.75),
        axis.line = element_line(colour = "black", size = 0.8))+
  labs(x = "", y = "Relative Abundance (%)", fill = "",size = "12", face = "bold")+
  scale_fill_manual(values = TRT_color)+
  scale_y_continuous(breaks=seq(0,10,1))+
  theme(legend.position = "none")+
  geom_text(aes(x=Sample, y= Value, label= round(Value, digits = 3)), hjust = -0.7,
            position=position_dodge(-0.8), size = 4.0,
            fontface = 2, family = "Times New Roman")
p1
ggsave("231010_ANCOM-BC_L6_Relative_abundance.png",
       width = 9, height = 9, units = "in", dpi = 600)

# ANCOM-BC Plot
# L2 and L6
## Input
setwd("D:/220729_Metagenomics_R_script/220808_R/230911_Hyunho_Soil_Metal/ANCOMBC2/231010_Plot/231010_L6/")
A <- read.table("231010_Plot_L2_L6_Ancombc_Holm_Bonf_Sig.txt", sep = "\t", header = TRUE)

# CON vs. MS
CON.MS <- A %>% subset(A$Variable == "CON vs. MS")
CON.MS$col <- gsub("blue","#3988C6",CON.MS$col)
CON.MS$col <- gsub("red1","#f65A46",CON.MS$col)
name=c("1 : plane","2 : bold","3 : italic","4 : bold and italic")
g = ggplot(CON.MS, aes(x = LFC, y = Taxa, fill = Variable, color = Variable))+
  geom_errorbarh(aes(xmin=Lower,xmax=Upper),color="black", height = 0.1) +
  geom_point(shape = 18, size = 6.0, color = CON.MS$col)+
  scale_x_continuous(breaks = seq(-5,6.5,1))+
  geom_vline(xintercept = 0, color="grey55", size = 1.2, linetype = "dashed")+
  ggtitle("")+
  theme_grey()
g
g1 <- g + theme(legend.title = element_blank(), 
                legend.text = element_blank(),
                legend.position = "none",
                axis.title.y = element_blank(),
                axis.text.x = element_text(face = "bold",colour = "black", size = 15, family = "Times New Roman"), 
                axis.text.y = element_text(face = "bold",colour = "black", size = 15, family = "Times New Roman"), 
                axis.title.x = element_text(face = "bold", size = 16, colour = "black", family = "Times New Roman",
                                            margin = margin(t = 0, r = 0, b = 25, l = 0)),
                strip.text.x = element_text(face = "bold", size = 15, colour = "black", family = "Times New Roman"),
                plot.title = element_text(face = "bold", size = 18, colour = "black", family = "Times New Roman"),
                axis.line = element_line(colour = "black", size = 0.8))+ 
  labs(x= "", y = "", fill = "")+
  coord_cartesian(xlim = c(-4.0,6.25), ylim = c(1,26), clip = "off")+ # change ylim number matching the number of taxa
  annotate("text", x = 3.5, y = -2.4, label = "Increased in MS",
           family = "Times New Roman", fontface = 2, size = 5.0)+
  annotate("text", x = 0.0, y = -1.1, label = "Log fold change",
           family = "Times New Roman", fontface = 2, size = 5.5)+
  annotate("segment", x = 0.5, y = -1.9, xend = 6.0, yend = -1.9, 
           size = 1.4, lineend = "round", linejoin = "round",
           arrow = arrow(length = unit(0.15, "inches")))+
  annotate("text", x = -2.0, y = -2.4, label = "Increased in CON",
           family = "Times New Roman", fontface = 2, size = 5.0)+
  annotate("segment", x = -0.5, y = -1.9, xend = -4.0, yend = -1.9, 
           size = 1.4, lineend = "round", linejoin = "round",
           arrow = arrow(length = unit(0.15, "inches")))+
  scale_y_discrete(limits = c("UCG_Rhodospirillales","g_Anaeromyxobacter","UCG_Cyclobacteriaceae","g_Afipia",
                              "g_B1-7BS","g_Pricia","g_Magnetospirillaceae","UCG_Elsterales","g_Herminiimonas",
                              "g_Ureibacillus","g_Nitrospira","UCG_Anaerolineaceae","g_AKAU4049","g_AKYG1722",
                              "g_Nannocystis","g_JG30-KF-CM66","UG_Gemmatimonadaceae","g_Nocardioides",
                              "g_Skermanella","UG_Sphingomonadaceae","g_Gitt-GS-136",
                              "g_Longimicrobiaceae","g_Luteimonas","UG_Micrococcaceae", "p_Hydrogenedentes","p_Bacteroidota"),
                   labels = c(expression(bold("UCG_Rhodospirillales")),
                              expression(bold("g_")*bolditalic("Anaeromyxobacter")),
                              expression(bold("UCG_Cyclobacteriaceae")),
                              expression(bold("g_")*bolditalic("Afipia")),
                              expression(bold("g_B1-7BS")),
                              expression(bold("g_")*bolditalic("Pricia")),
                              expression(bold("g_")*bolditalic("Magnetospirillaceae")),
                              expression(bold("UCG_Elsterales")),
                              expression(bold("g_")*bolditalic("Herminiimonas")),
                              expression(bold("g_")*bolditalic("Ureibacillus")),
                              expression(bold("g_")*bolditalic("Nitrospira")),
                              expression(bold("UCG_Anaerolineaceae")),
                              expression(bold("g_AKAU4049")),
                              expression(bold("g_AKYG1722")),
                              expression(bold("g_")*bolditalic("Nannocystis")),
                              expression(bold("g_JG30-KF-CM66")),
                              expression(bold("UG_Gemmatimonadaceae")),
                              expression(bold("g_")*bolditalic("Nocardioides")),
                              expression(bold("g_")*bolditalic("Skermanella")),
                              expression(bold("UG_Sphingomonadaceae")),
                              expression(bold("g_Gitt-GS-136")),
                              expression(bolditalic("Longimicrobiaceae")),
                              expression(bold("g_")*bolditalic("Luteimonas")),
                              expression(bold("UG_Micrococcaceae")),
                              expression(bold("p_Hydrogenedentes")),
                              expression(bold("p_Bacteroidota"))))
g1

p_right <- ggplot(data = CON.MS, aes(x=0, y=Taxa, label=P)) +
  theme_void()+
  geom_text(aes(x = 0, y = Taxa, label = P),size = 5.2,
            hjust = 0.5, family = "Times New Roman", fontface="bold") +
  ggtitle(expression(atop(paste(bold("Adjusted")~bolditalic("P")))))+
  theme(plot.title = element_text(hjust = 0.5, vjust = -6, size = 15, face = 2, family = "Times New Roman"),
        legend.text = element_blank(),
        legend.position = "none")
p_right

library("patchwork")
layout <- c(
  area(t = 0, l = 0, b = 10, r = 6.9), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  area(t = 0, l = 6.0, b = 10, r = 8) # right most plot starts at top of page, begins where middle plot ends (l=9, and middle plot is r=9), goes to bottom of page (b=30), and extends two units wide (r=11)
)
# final plot arrangement
P_CON.MS <- g1 + theme(plot.title = element_text(hjust = -2.1)) + p_right + plot_layout(design = layout)
P_CON.MS
ggsave("231010_ANCOM-BC_L2_L6_CON_MS.png",
       width = 10, height = 10, units = "in", dpi = 600)

# L6 Relative abundance bar plot
library("plyr")
library(reshape2)
B <-  read.table("231010_L2_L6_bar_plot_Relative_abundance.txt", sep="\t", header = T)
# Change the format of table for plotting
colnames(B) <- c("Sample","CON","MS")
B
AM = melt(B, id = c("Sample"))
AM
AM$variable <- factor(AM$variable, levels = c("CON","MS"))
df=ddply(AM, c("Sample","variable"), summarise, mean=mean(value), sd=sd(value), n= length(value), se=sd/sqrt(n))
df
colnames(df)[3]="Value"
df
mx <-  ggplot(df, aes(Sample, Value, fill=variable))+
  theme_gray()

# font --> face = 1 ---
par(mfrow=c(2,2))
name=c("1 : plane","2 : bold","3 : italic","4 : bold and italic")
TRT_color <- c("CON"="#f6d55c", "MS"="#80A0FF")
# phylum to genus
p1=mx+geom_bar(stat="identity", position=position_dodge(-0.8), width = 0.8)+
  coord_flip(ylim = c(0,12))+
  scale_x_discrete(limits = c("UCG_Rhodospirillales","g_Anaeromyxobacter","UCG_Cyclobacteriaceae","g_Afipia",
                              "g_B1-7BS","g_Pricia","g_Magnetospirillaceae","UCG_Elsterales","g_Herminiimonas",
                              "g_Ureibacillus","g_Nitrospira","UCG_Anaerolineaceae","g_AKAU4049","g_AKYG1722",
                              "g_Nannocystis","g_JG30-KF-CM66","UG_Gemmatimonadaceae","g_Nocardioides",
                              "g_Skermanella","UG_Sphingomonadaceae","g_Gitt-GS-136",
                              "g_Longimicrobiaceae","g_Luteimonas","UG_Micrococcaceae", "p_Hydrogenedentes","p_Bacteroidota"),
                   labels = c(expression(bold("UCG_Rhodospirillales")),
                              expression(bold("g_")*bolditalic("Anaeromyxobacter")),
                              expression(bold("UCG_Cyclobacteriaceae")),
                              expression(bold("g_")*bolditalic("Afipia")),
                              expression(bold("g_B1-7BS")),
                              expression(bold("g_")*bolditalic("Pricia")),
                              expression(bold("g_")*bolditalic("Magnetospirillaceae")),
                              expression(bold("UCG_Elsterales")),
                              expression(bold("g_")*bolditalic("Herminiimonas")),
                              expression(bold("g_")*bolditalic("Ureibacillus")),
                              expression(bold("g_")*bolditalic("Nitrospira")),
                              expression(bold("UCG_Anaerolineaceae")),
                              expression(bold("g_AKAU4049")),
                              expression(bold("g_AKYG1722")),
                              expression(bold("g_")*bolditalic("Nannocystis")),
                              expression(bold("g_JG30-KF-CM66")),
                              expression(bold("UG_Gemmatimonadaceae")),
                              expression(bold("g_")*bolditalic("Nocardioides")),
                              expression(bold("g_")*bolditalic("Skermanella")),
                              expression(bold("UG_Sphingomonadaceae")),
                              expression(bold("g_Gitt-GS-136")),
                              expression(bolditalic("Longimicrobiaceae")),
                              expression(bold("g_")*bolditalic("Luteimonas")),
                              expression(bold("UG_Micrococcaceae")),
                              expression(bold("p_Hydrogenedentes")),
                              expression(bold("p_Bacteroidota"))))+
  geom_errorbar(aes(ymin=Value-se, ymax=Value+se),
                position=position_dodge(-0.8), width=0.1)+
  ggtitle("")+
  theme(axis.text.x = element_text(angle = 0, size = 15, family = "Times New Roman", colour = "black", face = "bold"), 
        axis.title.x = element_text(size = 15, family = "Times New Roman", face = "bold"),
        axis.title.y = element_text(size = 15, family = "Times New Roman", face = "bold"), 
        legend.title = element_blank(),
        legend.text = element_blank(),
        axis.text.y = element_text(colour = "black", family = "Times New Roman", size = 15, face=2),
        plot.title = element_text(size = 18, family = "Times New Roman", face = "bold"),
        text = element_text(size = 15, family = "Times New Roman", face = "bold", colour = "black",
                            vjust = 0.75),
        axis.line = element_line(colour = "black", size = 0.8))+
  labs(x = "", y = "Relative Abundance (%)", fill = "",size = "12", face = "bold")+
  scale_fill_manual(values = TRT_color)+
  scale_y_continuous(breaks=seq(0,12,1))+
  theme(legend.position = "none")+
  geom_text(aes(x=Sample, y= Value, label= round(Value, digits = 3)), hjust = -0.7,
            position=position_dodge(-0.8), size = 4.0,
            fontface = 2, family = "Times New Roman")
p1
ggsave("231010_ANCOM-BC_L2_L6_Relative_abundance.png",
       width = 9, height = 9, units = "in", dpi = 600)