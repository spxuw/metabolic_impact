############################################################
# gut microbiome-plasma metabilites impact pipeline
# output:
#   1. Figure 7: acvd


library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)


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
colnames(diff_mat_null_test) <- colnames(diff_mat_test)

# Pick one representative metabolite for panel A
metab_focus <- c("dimethylglycine","N-acetylserine","N-acetylmethionine")
print(intersect(colnames(diff_mat_test),metab_focus))

# =========================
# 1) COMBINE DISEASE DATA
# =========================
meta_dis <- meta_test

# Ensure row alignment
stopifnot(nrow(meta_dis) == nrow(diff_mat_test))
stopifnot(nrow(meta_dis) == nrow(diff_mat_null_test))

# Build long table for selected metabolites
df_long <- bind_rows(lapply(metab_focus, function(met) {
  data.frame(
    genus = meta_dis$V4,
    sample = meta_dis$V5,
    metabolite = met,
    impact = diff_mat_test[, met],
    null = diff_mat_null_test[, met],
    stringsAsFactors = FALSE
  )
}))

# Absolute perturbation magnitude
df_long <- df_long %>%
  mutate(
    impact_abs = abs(impact),
    null_abs   = abs(null),
    gain       = impact_abs - null_abs
  )


# =========================
# 2) PANEL A
# Top genera for one representative metabolite
# =========================
df_genus_focus <- df_long %>%
  filter(metabolite == metab_focus) %>%
  group_by(genus) %>%
  summarise(
    mean_impact = mean(impact_abs, na.rm = TRUE),
    mean_null   = mean(null_abs, na.rm = TRUE),
    mean_gain   = mean(gain, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 10) %>%                       # optional stability filter
  arrange(desc(mean_gain))

top_genus_focus <- df_genus_focus %>%
  slice_head(n = 15)

df_plot_A <- top_genus_focus %>%
  select(genus, mean_impact, mean_null) %>%
  pivot_longer(
    cols = c(mean_impact, mean_null),
    names_to = "model",
    values_to = "value"
  ) %>%
  mutate(
    model = recode(model,
                   mean_impact = "Interaction-aware",
                   mean_null   = "Null")
  )

df_plot_A$genus <- factor(
  df_plot_A$genus,
  levels = top_genus_focus$genus[order(top_genus_focus$mean_gain)]
)

gA <- ggplot(df_plot_A, aes(x = genus, y = value, color = model, group = genus)) +
  geom_segment(
    data = top_genus_focus,
    aes(x = genus, xend = genus, y = mean_null, yend = mean_impact),
    inherit.aes = FALSE,
    color = "grey70"
  ) +
  geom_point(position = position_dodge(width = 0.25), size = 2) +
  custom_theme1 + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + 
  scale_color_manual(values = c("Interaction-aware" = "#d73027", "Null" = "#4575b4")) +
  xlab("") +
  ylab(paste0("Mean absolute perturbation\n(", metab_focus, ")")) +
  ggtitle("Representative ACVD-related metabolite")


# =========================
# 4) PANEL B
# Sample-level effect: interaction-aware vs null
# =========================
df_sample_sum <- df_long %>%
  group_by(sample, metabolite) %>%
  summarise(
    impact_abs = mean(impact_abs, na.rm = TRUE),
    null_abs   = mean(null_abs, na.rm = TRUE),
    gain       = mean(gain, na.rm = TRUE),
    .groups = "drop"
  )

gB <- ggplot(df_sample_sum, aes(x = null_abs, y = impact_abs)) +
  geom_point(alpha = 0.5, size = 1.5, color = "#80b1d3") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ metabolite, scales = "free") +
  custom_theme1 +
  xlab("Null perturbation magnitude") +
  ylab("Interaction-aware perturbation magnitude") +
  ggtitle("Sample-level metabolite responses")

# =========================
# 5) COMBINE
# =========================
p_final <- ggarrange(
  gA, gB, 
  ncol = 1, nrow = 2,
  heights = c(1.1, 1.0),
  labels = c("a", "b")
)

ggsave(
  p_final,
  file = "../figures/acvd_metabolite_tuning_interaction_vs_null.pdf",
  width = 8.5, height = 7.6, scale = 0.8
)

library(ComplexHeatmap)
library(circlize)

coef_focus <- coef[top_genus_focus$genus,metab_focus]





