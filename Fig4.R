############################################################
# gut microbiome-plasma metabilites impact pipeline
# output:
#   1. Figure 5: heatmap
#   2. Figure S3: prediction error


library(ggpubr)
library(ggridges)
library(reshape2)
library(ggrepel)
library(stringr)


setwd("/Users/xuwenwang/Library/CloudStorage/Dropbox/Projects/KeyPlasma/code")

rm(list = ls())

custom_theme <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(1.5, "mm"),
    axis.ticks = element_line(linewidth = 0.2, color = "black"),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.3, color = "black"),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    strip.text = element_text(face = "bold.italic", size = 8),
    strip.background = element_blank(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )

custom_theme1 <- theme_bw() +
  theme(
    # panel.grid.major = element_blank(),
    #  panel.grid.minor = element_blank(),
    axis.ticks.length = unit(1.5, "mm"),
    axis.ticks = element_line(linewidth = 0.2, color = "black"),
    #panel.border = element_blank(),
    #axis.line = element_line(linewidth = 0.3, color = "black"),
    #axis.line.x.top = element_blank(),
    #axis.line.y.right = element_blank(),
    strip.text = element_text(face = "bold.italic", size = 8),
    strip.background = element_blank(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )

M <- read.csv(file = "../results/KeyPlasma/qtst1.csv",header = F)
N <- read.csv(file = "../data/KeyPlasma/P_train.csv",header = F)
meta <- read.csv(file = "../data/KeyPlasma/meta_trainr.csv",header = F)
coef <- read.csv(file = "../data/KeyPlasma/assoc_mat.csv",header = T, row.names = 1, check.names = F)

M1 <- as.matrix(M) %*% as.matrix(coef)
N1 <- as.matrix(t(N)) %*% as.matrix(coef)

mu <- colMeans(N1, na.rm = TRUE)
sdv <- apply(N1, 2, sd, na.rm = TRUE)
sdv[sdv == 0] <- 1

M1 <- scale(M1, center = mu, scale = sdv)
N1 <- scale(N1, center = mu, scale = sdv)

print(dim(meta))
print(dim(M1))
print(dim(N1))

diff_mat <- t(vapply(seq_len(nrow(meta)), function(i) {
  m_col <- N1[meta$V3[i], ]
  n_row <- M1[i, ]
  n_row - m_col
}, numeric(ncol(M1))))

diff_mat_null <- t(vapply(seq_len(nrow(meta)), function(i) {
  m_col <- N[, meta$V3[i]]
  m_col[meta$V2[i]] <- 0
  n_row <- N[, meta$V3[i]]
  
  m_col1 <- as.matrix(t(m_col)) %*% as.matrix(coef)
  n_row1 <- as.matrix(t(n_row)) %*% as.matrix(coef)
  
  m_col1 <- (m_col1 - mu) / sdv
  n_row1 <- (n_row1 - mu) / sdv
  
  as.numeric(n_row1 - m_col1)
}, numeric(ncol(M1))))


M <- read.csv(file = "../results/KeyPlasma/qtst2.csv",header = F)
N <- read.csv(file = "../data/KeyPlasma/P_test.csv",header = F)
meta_test <- read.csv(file = "../data/KeyPlasma/meta_testr.csv",header = F)

M1 <- as.matrix(M) %*% as.matrix(coef)
N1 <- as.matrix(t(N)) %*% as.matrix(coef)

M1 <- scale(M1, center = mu, scale = sdv)
N1 <- scale(N1, center = mu, scale = sdv)

print(dim(meta_test))
print(dim(M1))
print(dim(N1))

diff_mat_test <- t(vapply(seq_len(nrow(meta_test)), function(i) {
  m_col <- N1[meta_test$V3[i], ]
  n_row <- M1[i, ]
  n_row - m_col
}, numeric(ncol(M1))))


diff_mat_null_test <- t(vapply(seq_len(nrow(meta_test)), function(i) {
  m_col <- N[, meta_test$V3[i]]
  m_col[meta_test$V2[i]] <- 0
  n_row <- N[, meta_test$V3[i]]
  
  m_col1 <- as.matrix(t(m_col)) %*% as.matrix(coef)
  n_row1 <- as.matrix(t(n_row)) %*% as.matrix(coef)
  
  m_col1 <- (m_col1 - mu) / sdv
  n_row1 <- (n_row1 - mu) / sdv
  
  as.numeric(n_row1 - m_col1)
}, numeric(ncol(M1))))

gain_mat <- abs(rbind(diff_mat,diff_mat_test)) - abs(rbind(diff_mat_null,diff_mat_null_test))


disease_label <- read.csv(file = "../data/KeyPlasma/stool_wgs_sample_table.csv")
disease_label$V5 <- rownames(disease_label)

meta_combined <- disease_label
meta_combined$disease_main <- sapply(strsplit(as.character(meta_combined$disease), ";"), `[`, 1)

meta_removal <- rbind(meta,meta_test)


meta_removal2 <- meta_removal %>%
  dplyr::left_join(disease_label, by = "V5") %>%
  mutate(disease_main = sapply(strsplit(as.character(disease), ";"), `[`, 1))


sample_ids <- meta_removal2$V5
idx_by_sample <- split(seq_len(nrow(gain_mat)), sample_ids)

sample_metab_mat <- t(sapply(idx_by_sample, function(idx) {
  colMeans(gain_mat[idx, , drop = FALSE], na.rm = TRUE)
}))
sample_metab_mat <- as.matrix(sample_metab_mat)


meta_combined <- meta_combined[match(rownames(sample_metab_mat), rownames(meta_combined)), ]
df_sample_metab <- cbind(meta_combined, as.data.frame(sample_metab_mat))

keep_dis <- c("healthy", "ACVD", "T2D")

meta_sample_sub <- meta_combined %>%
  filter(disease_main %in% keep_dis)

sample_metab_sub <- sample_metab_mat[rownames(meta_sample_sub), , drop = FALSE]

mean_by_disease <- function(dis) {
  idx <- meta_sample_sub$disease_main == dis
  colMeans(sample_metab_sub[idx, , drop = FALSE], na.rm = TRUE)
}

gain_healthy <- mean_by_disease("healthy")
gain_acvd    <- mean_by_disease("ACVD")
gain_t2d     <- mean_by_disease("T2D")

df_metab_wide <- data.frame(
  metabolite = colnames(sample_metab_sub),
  healthy = gain_healthy,
  ACVD = gain_acvd,
  T2D = gain_t2d
) %>%
  mutate(
    gain_acvd = ACVD - healthy,
    gain_t2d  = T2D - healthy
  )

top_n <- 20

df_top <- bind_rows(
  df_metab_wide %>%
    arrange(desc(gain_acvd)) %>%
    slice_head(n = top_n) %>%
    transmute(metabolite, contrast = "ACVD vs healthy", gain = gain_acvd),
  
  df_metab_wide %>%
    arrange(desc(gain_t2d)) %>%
    slice_head(n = top_n) %>%
    transmute(metabolite, contrast = "T2D vs healthy", gain = gain_t2d)
)

g1 <- ggplot(df_top[df_top$contrast=="ACVD vs healthy",], aes(x = reorder(metabolite, gain), y = gain, fill = contrast)) +
  geom_col(show.legend = FALSE,fill="#7BAFD9") +
  coord_flip() +
  custom_theme1 +
  xlab("") +
  ylab("Interaction-specific metabolite perturbation")

g2 <- ggplot(df_top[df_top$contrast=="T2D vs healthy",], aes(x = reorder(metabolite, gain), y = gain, fill = contrast)) +
  geom_col(show.legend = FALSE,fill="#F8968C") +
  coord_flip() +
  custom_theme1 +
  xlab("") +
  ylab("Interaction-specific metabolite perturbation")


annot <- read.csv(file = "/Users/xuwenwang/Library/CloudStorage/Dropbox/Projects/KeyPlasma/data/Supplementary Table 2.tsv", sep = "\t")
df_metab_wide$subclass <- annot$Metabolite.subclass[match(df_metab_wide$metabolite,annot$Metabolite)]
df_metab_wide$class <- annot$Metabolite.class[match(df_metab_wide$metabolite,annot$Metabolite)]


top_n <- 20

top_metabs <- df_metab_wide %>%
  arrange(desc(abs(gain_acvd))) %>%
  slice_head(n = top_n) %>%
  pull(metabolite)

top_metabs2 <- df_metab_wide %>%
  arrange(desc(abs(gain_t2d))) %>%
  slice_head(n = top_n) %>%
  pull(metabolite)


top_metabs <- unique(c(top_metabs1,top_metabs2))


cor_mat <- cor(sample_metab_mat, method = "spearman")
cor_sub <- cor_mat[top_metabs, top_metabs]

subclass_vec <- df_metab_wide$subclass[
  match(top_metabs, df_metab_wide$metabolite)
]
names(subclass_vec) <- top_metabs

class_vec <- df_metab_wide$class[
  match(top_metabs, df_metab_wide$metabolite)
]
names(class_vec) <- top_metabs

library(ComplexHeatmap)
library(circlize)

subclass_levels <- unique(subclass_vec)

col_subclass <- structure(
  colorRampPalette(c("#4DBBD5", "#E64B35", "#00A087", "#3C5488", "#F39B7F"))(length(subclass_levels)),
  names = subclass_levels
)

gain_acvd_vec <- df_metab_wide$gain_acvd[
  match(top_metabs, df_metab_wide$metabolite)
]

gain_t2d_vec <- df_metab_wide$gain_t2d[
  match(top_metabs, df_metab_wide$metabolite)
]

names(gain_acvd_vec) <- top_metabs
names(gain_t2d_vec)  <- top_metabs

col_gain <- colorRamp2(
  c(min(gain_acvd_vec, gain_t2d_vec, na.rm = TRUE), 0,
    max(gain_acvd_vec, gain_t2d_vec, na.rm = TRUE)),
  c("blue", "white", "red")
)

ha_row <- rowAnnotation(
  subclass = subclass_vec,
  ACVD_gain = gain_acvd_vec,
  T2D_gain  = gain_t2d_vec,
  
  col = list(
    subclass = col_subclass,
    ACVD_gain = col_gain,
    T2D_gain  = col_gain
  ),
  
  annotation_name_gp = gpar(fontsize = 8)
)

pdf("../figures/acvd_t2d_heatmap.pdf",width = 12.8,height = 8)
Heatmap(
  cor_sub,
  name = "Spearman rho",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  left_annotation = ha_row
)
dev.off()


## class level ...........................
class_levels <- unique(class_vec)

col_class <- structure(
  colorRampPalette(c("#4DBBD5", "#E64B35", "#00A087", "#3C5488", "#F39B7F"))(length(class_levels)),
  names = class_levels
)

ha_row <- rowAnnotation(
  subclass = class_vec,
  ACVD_gain = gain_acvd_vec,
  T2D_gain  = gain_t2d_vec,
  
  col = list(
    subclass = col_class,
    ACVD_gain = col_gain,
    T2D_gain  = col_gain
  ),
  
  annotation_name_gp = gpar(fontsize = 8)
)

pdf("../figures/acvd_t2d_heatmap_class.pdf",width = 10,height = 8)
Heatmap(
  cor_sub,
  name = "Spearman rho",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  left_annotation = ha_row
)
dev.off()


# A shared metabolite module exhibits coordinated behavior but is differentially regulated by ecological interactions across diseases, 
# highlighting disease-specific rewiring of microbiome–metabolome coupling.

# compare the inconsitance of null and interaction arare

qtrn <- read.csv(file = "../results/KeyPlasma/qtrn.csv",header = F)
ptrn <- read.csv(file = "../data/KeyPlasma/P_train.csv",header = F)
ptrn <- as.data.frame(t(ptrn))

bc_dis <- c()
for (i in 1:nrow(qtrn)){
  print(i)
  bc_dis <- c(bc_dis, cor(as.numeric(ptrn[i,]),as.numeric(qtrn[i,])))
}

med_val <- median(bc_dis, na.rm = TRUE)

gc <- ggplot(data = data.frame(cor = bc_dis), aes(cor)) + 
  geom_histogram(color="white",fill="#2b8cbe") + 
  geom_vline(xintercept = med_val, linetype = "dashed", size = 1) +
  xlab("Pearson correlation between true and prediction") + 
  annotate("text",
           x = med_val,
           y = Inf,
           label = paste0("Median = ", round(med_val, 3)),
           vjust = 1.5,
           hjust = 1.1) +
  custom_theme

ggsave(gc,file="../figures/prediction_cnode.pdf",width=5, height=5,scale = 0.8)



