############################################################
# gut microbiome-plasma metabilites impact pipeline
# output:
#   1. Figure 8: drug intervention


library(ggpubr)
library(ggridges)
library(reshape2)
library(ggrepel)
library(stringr)
library(dplyr)
library(tidyr)
library(glmnet)
library(pROC)
library(purrr)


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


disease_label <- read.csv("../data/KeyPlasma/stool_wgs_sample_table.csv")
disease_label$sample_id <- rownames(disease_label)

disease_label$disease_main <- sapply(strsplit(as.character(disease_label$disease), ";"), `[`, 1)


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


library(readxl)
library(matrixStats)
library(reshape2)


otu_fiber <- read.csv(file = "../data/art_variation_public/data/data_fiber/otu_table_1000556.txt",skip = 1, sep = "\t")
meta_fiber <- read.csv(file = "../data/art_variation_public/data/data_fiber/metadata_29307330.txt",sep = "\t")

otu_genus <- otu_fiber %>%
  mutate(
    Genus = str_extract(taxonomy, "g__[^;]+") %>%
      str_replace("g__", "") %>%
      replace_na("Unknown")
  ) %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

otu_genus <- otu_genus %>%
  mutate(across(-Genus, ~ .x / sum(.x, na.rm = TRUE)))

otu_long <- otu_genus %>%
  pivot_longer(
    cols = -Genus,
    names_to = "Run",
    values_to = "abundance"
  )

df <- otu_long %>%
  left_join(meta_fiber, by = "Run")

df <- df %>%
  mutate(
    time = ifelse(SAMPLE_TYPE == "Post-intervention", "post", "pre")
  )

df_delta <- df %>%
  group_by(Subject, intervention, Genus, time) %>%
  summarise(abundance = mean(abundance), .groups = "drop") %>%
  pivot_wider(names_from = time, values_from = abundance) %>%
  mutate(delta = post - pre)


df_diff <- df_delta %>%
  group_by(Genus, intervention) %>%
  summarise(mean_delta = mean(delta, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = intervention, values_from = mean_delta) %>%
  mutate(
    delta_true = Prebiotic - Placebo
  )


# attach genus
impact_df <- as.data.frame(diff_mat)
impact_df$genus <- meta$V4

impact_null_df <- as.data.frame(diff_mat_null)
impact_null_df$genus <- meta$V4

# aggregate
impact_genus <- impact_df %>%
  group_by(genus) %>%
  summarise(across(everything(), median, na.rm = TRUE))

impact_null_genus <- impact_null_df %>%
  group_by(genus) %>%
  summarise(across(everything(), median, na.rm = TRUE))

impact_mat <- as.matrix(impact_genus[,-1])
rownames(impact_mat) <- impact_genus$genus

impact_null <- as.matrix(impact_null_genus[,-1])
rownames(impact_null) <- impact_null_genus$genus
colnames(impact_null) <- colnames(impact_mat)

common_genus <- intersect(df_diff$Genus,rownames(impact_mat))
common_pos <- intersect(common_genus,df_diff$Genus[df_diff$delta_true>0])
common_neg <- intersect(common_genus,df_diff$Genus[df_diff$delta_true<0])

impact_mat_common <- impact_mat[common_pos,colSums(abs(impact_mat))>150]

pdf("../figures/prebitoc_pos.pdf",width = 9.5,height = 5)
pheatmap(t(impact_mat_common),scale="column",fontsize = 8)
dev.off()

impact_mat_common <- impact_mat[common_neg,colSums(abs(impact_mat))>150]
pdf("../figures/prebitoc_neg.pdf",width = 9.5,height = 5)
pheatmap(t(impact_mat_common),scale="column",fontsize = 8)
dev.off()


# mean effect per metabolite
met_pos <- colMeans(impact_mat[common_pos, ], na.rm = TRUE)
met_neg <- colMeans(impact_mat[common_neg, ], na.rm = TRUE)

df_imp <- data.frame(
  metabolite = names(met_pos),
  diff = met_pos - met_neg
) %>%
  mutate(abs_diff = abs(diff)) %>%
  arrange(desc(abs_diff))

top10_imp <- df_imp %>% slice_head(n = 10)


# null model
met_pos_null <- colMeans(impact_null[common_pos, ], na.rm = TRUE)
met_neg_null <- colMeans(impact_null[common_neg, ], na.rm = TRUE)

df_null <- data.frame(
  metabolite = names(met_pos_null),
  diff = met_pos_null - met_neg_null
) %>%
  mutate(abs_diff = abs(diff)) %>%
  arrange(desc(abs_diff))

top10_null <- df_null %>% slice_head(n = 10)

df_plot <- bind_rows(
  top10_imp %>% mutate(model = "Interaction"),
  top10_null %>% mutate(model = "Null")
)

df_plot$metabolite <- factor(df_plot$metabolite,
                             levels = unique(df_plot$metabolite))

g_top <- ggplot(df_plot,
                aes(x = diff,
                    y = metabolite,
                    color = diff>0)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~model, scales = "free_y") +
  custom_theme1 +
  xlab("Difference (POS - NEG impact)") +
  ylab("Top discriminative metabolites")+
  theme(legend.position = "none")


ggsave(g_top,file="../figures/prebiotic_top.pdf",width=9, height=5,scale = 0.8)

