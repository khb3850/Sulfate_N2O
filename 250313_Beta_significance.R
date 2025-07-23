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
setwd("/Users/khb3850/Library/CloudStorage/GoogleDrive-khb3850@gmail.com/My Drive/Code/PostDoc/250305_Hyunho_PCA/250313_Soil_Beta_diversity_Update")

# Data importing from qza result file of beta diversity (pcoa_results.qza)
metadata<-read_q2metadata("230911_Meta_data.txt")
Unwe_uni_frac<-read_qza("unweighted_unifrac_pcoa_results.qza")
we_uni_frac<-read_qza("weighted_unifrac_pcoa_results.qza")

# PCoA 좌표 추출
unwe_pcoa_df <- Unwe_uni_frac$data$Vectors %>% select(SampleID, PC1, PC2, PC3)
we_pcoa_df <- we_uni_frac$data$Vectors %>% select(SampleID, PC1, PC2, PC3)

# 메타데이터와 병합
unwe_pcoa_df <- left_join(unwe_pcoa_df, metadata, by = c("SampleID" = "SampleID"))
we_pcoa_df <- left_join(we_pcoa_df, metadata, by = c("SampleID" = "SampleID"))

# 📥 거리 행렬 가져오기
unwe_dist <- read_qza("250331_Beta_Soil/unweighted_unifrac_distance_matrix.qza")$data
we_dist <- read_qza("250331_Beta_Soil/weighted_unifrac_distance_matrix.qza")$data
unwe_dist <- as.dist(unwe_dist)
we_dist <- as.dist(we_dist)

# PERMANOVA
## Unweighted UniFrac
set.seed(123)
adonis2(unwe_dist ~ TRT, data = metadata, permutations = 9999)
# Pairwise PERMANOVA
set.seed(123)
pairwise.adonis2(unwe_dist ~ TRT, data = metadata, nperm = 9999, p.adjust.methods = "bonferroni")

## Weighted UniFrac
set.seed(123)
adonis2(we_dist ~ TRT, data = metadata, permutations = 9999)
# Pairwise PERMANOVA
set.seed(123)
pairwise.adonis2(we_dist ~ TRT, data = metadata, nperm = 9999, p.adjust.methods = "bonferroni")

# PERMDISP
## Unweighted UniFrac
unwe_dispersion <- betadisper(unwe_dist, metadata$TRT)
permutest(unwe_dispersion, permutations = 9999)

## Weighted UniFrac
we_dispersion <- betadisper(we_dist, metadata$TRT)
permutest(we_dispersion, permutations = 9999)


# Bar-plot for PERMDISP
## unweighted UniFrac 거리 기준
disp <- betadisper(unwe_dist, metadata$TRT)
distance_df <- data.frame(
  SampleID = names(disp$distances),
  Group = metadata$TRT[match(names(disp$distances), metadata$SampleID)],
  DistanceToCentroid = disp$distances
)
color <- c("CON" = "#f6d55c", "MS" = "#FFA080", "ZVM" = "#80A0FF")

# 그룹 순서 고정 (중요!)
distance_df$Group <- factor(distance_df$Group, levels = c("CON", "MS", "ZVM"))

# 그림
ggplot(distance_df, aes(x = Group, y = DistanceToCentroid, fill = Group)) +
  geom_boxplot(alpha = 0.8, width = 0.6, outlier.shape = 21, outlier.fill = "black") +
  scale_fill_manual(values = color) +
  theme_grey(base_family = "Times New Roman") +
  ylim(0.25,0.4)+
  labs(
    title = "",
    y = "Distance to centroid", x = NULL
  ) +
  #annotate("text", x = 1, y = max(distance_df$DistanceToCentroid) + 0.05,
  #         label = "PERMDISP  F = 18.167\np < 0.001", hjust = 0, size = 5) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 20, face = "bold")
  )
ggsave("250403_Unweighted_PERMDISP.png", width = 6, height = 5.5, units = "in", dpi = 600)


## weighted UniFrac 거리 기준
disp <- betadisper(we_dist, metadata$TRT)
distance_df <- data.frame(
  SampleID = names(disp$distances),
  Group = metadata$TRT[match(names(disp$distances), metadata$SampleID)],
  DistanceToCentroid = disp$distances
)
color <- c("CON" = "#f6d55c", "MS" = "#FFA080", "ZVM" = "#80A0FF")

# 그룹 순서 고정 (중요!)
distance_df$Group <- factor(distance_df$Group, levels = c("CON", "MS", "ZVM"))

# 그림
ggplot(distance_df, aes(x = Group, y = DistanceToCentroid, fill = Group)) +
  geom_boxplot(alpha = 0.8, width = 0.6, outlier.shape = 21, outlier.fill = "black") +
  scale_fill_manual(values = color) +
  theme_grey(base_family = "Times New Roman") +
  ylim(0.0,0.10)+
  labs(
    title = "",
    y = "Distance to centroid", x = NULL
  ) +
  #annotate("text", x = 1, y = max(distance_df$DistanceToCentroid) + 0.05,
  #         label = "PERMDISP  F = 18.167\np < 0.001", hjust = 0, size = 5) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 20, face = "bold")
  )
ggsave("250403_weighted_PERMDISP.png", width = 6, height = 5.5, units = "in", dpi = 600)







## Pairwise comparison --> Results are something different with QIIME2
qiime2_style_permdisp <- function(dist_mat, metadata, group_column, permutations = 9999) {
  groups <- unique(metadata[[group_column]])
  combs <- combn(groups, 2, simplify = FALSE)
  
  results <- data.frame(Group1 = character(),
                        Group2 = character(),
                        F = numeric(),
                        p_value = numeric(),
                        stringsAsFactors = FALSE)
  
  for (pair in combs) {
    # 1. 두 그룹 샘플만 선택
    idx <- metadata[[group_column]] %in% pair
    sub_metadata <- metadata[idx, ]
    sample_ids <- sub_metadata$SampleID
    group_labels <- droplevels(factor(sub_metadata[[group_column]]))
    
    # 2. 전체 거리행렬에서 해당 샘플만 추출 (거리 고정)
    sub_dist <- as.matrix(dist_mat)[sample_ids, sample_ids]
    dist_obj <- as.dist(sub_dist)
    
    # 3. 관측된 F
    disp <- betadisper(dist_obj, group_labels)
    f_obs <- anova(disp)$"F value"[1]
    
    # 4. permutation (label만 shuffle)
    f_null <- numeric(permutations)
    set.seed(42)  # 재현성
    for (i in 1:permutations) {
      shuffled <- sample(group_labels)
      disp_perm <- betadisper(dist_obj, shuffled)
      f_null[i] <- anova(disp_perm)$"F value"[1]
    }
    
    # 5. empirical p-value 계산 (QIIME2 방식)
    p_val <- (sum(f_null >= f_obs) + 1) / (permutations + 1)
    
    results <- rbind(results, data.frame(
      Group1 = pair[1],
      Group2 = pair[2],
      F = f_obs,
      p_value = p_val
    ))
  }
  
  # BH (FDR) 보정
  results$p_adj <- p.adjust(results$p_value, method = "BH")
  return(results)
}

# dist_mat: distance matrix in dist or matrix form
# metadata: dataframe with SampleID and TRT columns

# ensure metadata is in correct order
metadata <- metadata[match(rownames(as.matrix(we_dist)), metadata$SampleID), ]

qiime2_pairwise <- qiime2_style_permdisp(
  dist_mat = unwe_dist,       # already a dist object
  metadata = metadata,
  group_column = "TRT",
  permutations = 9999
)

qiime2_pairwise



