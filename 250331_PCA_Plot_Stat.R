# PCA analysis with PERMANOVA
## https://rpubs.com/collnell/manova
## https://github.com/pmartinezarbizu/pairwiseAdonis
install.packages('vegan')
install.packages("tibble")
install.packages("ggplot2")
install.packages("ggfortify")
install.packages('devtools')
devtools::install_github("cmartin/ggConvexHull")
library(devtools)
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
install.packages("pairwiseAdonis")
install.packages("tibble")
library("pairwiseAdonis")
library(ggplot2)
library(ggfortify)
library(ggConvexHull)
library(vegan) ##Community ecology: ordination, disversity & dissimilarities
library(dplyr)
library(reshape2)
library(tibble)
install.packages("ggrepel")  # 설치가 안 되어 있다면
library(ggrepel)

# Making dataframe for PCA and Statistical analysis
## Set the working directory
getwd()
setwd("/Users/khb3850/Library/Mobile Documents/com~apple~CloudDocs/1_All/Code/PostDoc/250305_Hyunho_PCA/")

# Load metadata
metadata <- read.table("250305_Metadata.txt", sep = "\t", header = TRUE)

# Load KO table (row = features, col = samples)
library(tibble)
KO <- read.table("250305_Raw_data.txt", sep = "\t", header = TRUE, row.names = 1)
KO_t <- as.data.frame(t(KO))
KO_t <- rownames_to_column(KO_t, var = "SampleID")

# Merge with metadata
KOs <- merge(metadata, KO_t, by = "SampleID")
KOs[, -(1:4)] <- sapply(KOs[, -(1:4)], as.numeric)

# PERMANOVA
KOs.matrix <- as.matrix(KOs[, -(1:4)])
KOs.prop <- decostand(KOs.matrix, method = "total")
KOs.dist <- vegdist(KOs.prop, method = "bray")

set.seed(123)
adonis2(KOs.dist ~ TRT_2, data = KOs, permutations = 9999)

# Pairwise PERMANOVA
set.seed(123)
pairwise.adonis2(KOs.dist ~ TRT_2, data = KOs, nperm = 9999, p.adjust.methods = "bonferroni")

# PCA
## https://stats.stackexchange.com/questions/53/pca-on-correlation-or-covariance
## There are two types of standard for PCA analysis, Covariance matrix (scale = FALSE) and correlation matrix (scale = TRUE)
## For the normalized data, we have to use result from Covariance matrix
PCA_result <- prcomp(KOs[, -(1:4)], center = TRUE, scale. = FALSE)
PCA_result
# Print proportioin of variance / Check elbow point
plot(PCA_result, type = "l")
## Check the total proportion from PC1 to PC nubmer at elbow point
## General standard for elbow point = 90%
summary(PCA_result)

PCA_data <- as.data.frame(PCA_result$x)
PCA_data <- cbind(PCA_data, KOs[, 1:4])
PCA_data$TRT_2 <- as.factor(PCA_data$TRT_2)
PCA_data$TRT_3 <- as.factor(PCA_data$TRT_3)
PCA_data

# Making your Figure
p <- ggplot(PCA_data, aes(x = PC1, y = PC2, color = TRT_2, shape = TRT_3)) +
  geom_point(size = 5, alpha = 0.95) +
  stat_ellipse(aes(group = TRT_2, fill = TRT_2), geom = "polygon", level = 0.90, alpha = 0.05, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", size = 0.8) +
  scale_shape_manual(values = c(18, 17, 15)) +
  scale_color_manual(values = c("Cu" = "#F8766D", "Fe" = "#7CAE00", "Mn" = "#00BFC4", "Mo" = "#C77CFF", "Zn" = "blue", "CON" = "black")) + # 금속별 색상
  theme_grey() +
  theme(
    axis.text = element_text(size = 16, family = "Times New Roman", face = "bold"),
    axis.title = element_text(size = 16, face = "bold", family = "Times New Roman"),
    legend.title = element_text(size = 14, family = "Times New Roman", face = "bold"),
    legend.text = element_text(size = 14, family = "Times New Roman", face = "bold"),
    legend.position = "bottom"
  ) +
  guides(
    color = guide_legend(title = "Micronutrient types", nrow = 1),
    shape = guide_legend(
      title = "Application rates (mg kg⁻¹)",
      nrow = 1,
      title.theme = element_text(
        size = 14,
        family = "Times New Roman",
        face = "bold"
      )
    )
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical"  # 각 legend를 줄 바꿈해서 쌓아줌
  ) +
  labs(title = "", 
       x = "PC1 (54.5%)",
       y = "PC2 (25.4%)")
p
p1 <- p + coord_fixed(ratio = 1, xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0))
p1


# 1. PCA 결과에서 loadings 추출
loadings <- as.data.frame(PCA_result$rotation[, 1:2])  # PC1, PC2에 대한 loading
loadings$Feature <- rownames(loadings)
# Feature 이름 순서에 따라 label 재매핑
loadings$label <- dplyr::case_when(
  loadings$Feature == "NH4"  ~ "bold(NH[4]^'+')",
  loadings$Feature == "NO3"  ~ "bold(NO[3]^'-')",
  loadings$Feature == "N2O"  ~ "bold(N[2]*O)",
  loadings$Feature == "pH"   ~ "bold('Soil pH')",
  loadings$Feature == "amoA" ~ "bold(amoA)",
  loadings$Feature == "nirS" ~ "bold(nirS)",
  loadings$Feature == "nosZ" ~ "bold(nosZ)"
)

# 2. 화살표 스케일링 (너무 길거나 짧지 않도록)
scale_factor <- 1  # 필요시 조절
loadings$PC1 <- loadings$PC1 * scale_factor
loadings$PC2 <- loadings$PC2 * scale_factor

# 3. 기존 PCA plot에 geom_segment로 화살표 추가
p1 <- p + 
  coord_fixed(ratio = 1, xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0)) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray30", alpha = 0.6,
               inherit.aes = FALSE) +
  ggrepel::geom_text_repel(
    data = loadings,
    aes(x = PC1, y = PC2, label = label),
    parse = TRUE,
    size = 4,
    family = "Times New Roman",
    fontface = "bold",
    color = "gray30",
    inherit.aes = FALSE
  )

# 결과 출력
p1
ggsave("250713_PCA_KOs_Type_50-100.png", width = 9, height = 9, units = "in", dpi = 600)

p1 <- p + coord_fixed(ratio = 1, xlim = c(-10, 10), ylim = c(-10, 10))
p1
ggsave("250405_PCA_KOs_Type_50-400_2.png", width = 9, height = 9, units = "in", dpi = 600)
