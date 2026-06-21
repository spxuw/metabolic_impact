############################################################
# gut microbiome-plasma metabilites impact pipeline
# output:
#   1. Figure 2: microbe-metabolite association network distribution
#   2. Figure S2: nestness using different threshold
#   3. Figure S3: Distance


library(ggpubr)
library(ggridges)
library(reshape2)
library(bipartite)

setwd("/Users/xuwenwang/Library/CloudStorage/Dropbox/Projects/KeyPlasma/code")

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


coef <- read.csv(file = "../data/KeyPlasma/assoc_mat.csv",header = T, row.names = 1, check.names = F)
mat_bin <- as.matrix(coef)
mat_bin[mat_bin!=0] <- 1

row_deg <- rowSums(mat_bin)
col_deg <- colSums(mat_bin)

df_deg <- rbind(
  data.frame(degree = row_deg, side = "row"),
  data.frame(degree = col_deg, side = "column")
)

x <- df_deg$degree[df_deg$side == "column"]
lambda_hat <- mean(x)

binwidth <- 5

# define bins
breaks <- seq(min(x), max(x) + binwidth, by = binwidth)

# compute expected counts per bin
bin_centers <- head(breaks, -1) + binwidth/2

expected_counts <- sapply(1:(length(breaks)-1), function(i) {
  k_vals <- breaks[i]:(breaks[i+1]-1)
  sum(dpois(k_vals, lambda_hat)) * length(x)
})

df_fit <- data.frame(
  x = bin_centers,
  y = expected_counts
)

g5 <- ggplot(data.frame(degree = x), aes(x = degree)) +
  geom_histogram(binwidth = binwidth, fill="#5495CF", color="white") +
  xlab("Degree") + ylab("Count") +
  ggtitle(paste0("Poisson fit (lambda = ", round(lambda_hat, 2), ")")) +
  custom_theme


fit_nb <- fitdistr(x, "Negative Binomial")
size <- fit_nb$estimate["size"]
mu <- fit_nb$estimate["mu"]

expected_counts_nb <- sapply(1:(length(breaks)-1), function(i) {
  k_vals <- breaks[i]:(breaks[i+1]-1)
  sum(dnbinom(k_vals, size = size, mu = mu)) * length(x)
})

df_fit_nb <- data.frame(
  x = bin_centers,
  y = expected_counts_nb
)

g5 <- g5 +
  geom_line(data = df_fit_nb, aes(x = x, y = y),
            color = "darkgreen", linewidth = 1.2)



g6 <- ggplot(df_deg[df_deg$side=="row",], aes(x = degree)) +
  geom_histogram(bins = 30, fill="#9C3F5D",color="white") +
  facet_wrap(~side, scales = "free") +xlab("Degree") + ylab("Count")+
  custom_theme


pdf("../figures/incident_matrix.pdf",width = 8.5,height = 4)
Heatmap(
  mat_bin[rowSums(mat_bin)!=0,],
  col = colorRamp2(
    c(min(mat_bin), 0, max(mat_bin)),
    c("white", "white", "#ca0020")
  ),
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  use_raster = FALSE
)
dev.off()

print(nested(mat_bin[rowSums(mat_bin)!=0,], method = "NODF"))
# 38.41459
row_norms <- sqrt(rowSums(coef^2))
row_norms[row_norms == 0] <- 1
coef_row_norm <- coef / row_norms
row_dist <- as.matrix(dist(coef_row_norm, method = "euclidean"))

col_norms <- sqrt(colSums(coef^2))
col_norms[col_norms == 0] <- 1
coef_col_norm <- t(t(coef) / col_norms)
col_dist <- as.matrix(dist(t(coef_col_norm), method = "euclidean"))
row_dist <- reshape2::melt(row_dist)
col_dist <- reshape2::melt(col_dist)
df_l2 <- rbind(
  data.frame(l2 = row_dist$value, side = "row"),
  data.frame(l2 = col_dist$value, side = "column")
)

g7 <- ggplot(df_l2[df_l2$side=="column",], aes(x = l2)) +
  geom_histogram(bins = 30, fill="#5495CF",color="white") +
  facet_wrap(~side, scales = "free") +xlab("L2 distance") + ylab("Count")+
  custom_theme


g8 <- ggplot(df_l2[df_l2$side=="row",], aes(x = l2)) +
  geom_histogram(bins = 30, fill="#9C3F5D",color="white") +
  facet_wrap(~side, scales = "free") +xlab("L2 distance") + ylab("Count")+
  custom_theme

p2 <- ggarrange(g5,g6,g7,g8,nrow = 2,ncol = 2,labels = c("a","b","c","d"))
ggsave(p2,file="../figures/disttibution.pdf",width=8, height=8,scale = 0.7)


## ---------------------------------------------------------------------------

coef <- read.csv(file = "../data/KeyPlasma/assoc_mat.csv",header = T, row.names = 1, check.names = F)
mat_bin <- as.matrix(coef)
mat_bin[abs(mat_bin)<3] <- 0
mat_bin[mat_bin!=0] <- 1

row_deg <- rowSums(mat_bin)
col_deg <- colSums(mat_bin)

df_deg <- rbind(
  data.frame(degree = row_deg, side = "row"),
  data.frame(degree = col_deg, side = "column")
)

g5 <- ggplot(df_deg[df_deg$side=="column",], aes(x = degree)) +
  geom_histogram(bins = 30, fill="#5495CF",color="white") +
  facet_wrap(~side, scales = "free") + xlab("Degree") + ylab("Count")+
  custom_theme


g6 <- ggplot(df_deg[df_deg$side=="row",], aes(x = degree)) +
  geom_histogram(bins = 30, fill="#9C3F5D",color="white") +
  facet_wrap(~side, scales = "free") +xlab("Degree") + ylab("Count")+
  custom_theme


pdf("../figures/incident_matrix_3.pdf",width = 8.5,height = 4)
Heatmap(
  mat_bin[rowSums(mat_bin)!=0,],
  col = colorRamp2(
    c(min(mat_bin), 0, max(mat_bin)),
    c("white", "white", "#ca0020")
  ),
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  use_raster = FALSE
)
dev.off()

print(nested(mat_bin[rowSums(mat_bin)!=0,], method = "NODF"))
# 38.41459

row_norms <- sqrt(rowSums(coef^2))
row_norms[row_norms == 0] <- 1
coef_row_norm <- coef / row_norms
row_dist <- as.matrix(dist(coef_row_norm, method = "euclidean"))

col_norms <- sqrt(colSums(coef^2))
col_norms[col_norms == 0] <- 1
coef_col_norm <- t(t(coef) / col_norms)
col_dist <- as.matrix(dist(t(coef_col_norm), method = "euclidean"))
row_dist <- reshape2::melt(row_dist)
col_dist <- reshape2::melt(col_dist)

df_l2 <- rbind(
  data.frame(l2 = row_dist$value, side = "row"),
  data.frame(l2 = col_dist$value, side = "column")
)

g7 <- ggplot(df_l2[df_l2$side=="column",], aes(x = l2)) +
  geom_histogram(bins = 30, fill="#5495CF",color="white") +
  facet_wrap(~side, scales = "free") +xlab("L2 distance") + ylab("Count")+
  custom_theme


g8 <- ggplot(df_l2[df_l2$side=="row",], aes(x = l2)) +
  geom_histogram(bins = 30, fill="#9C3F5D",color="white") +
  facet_wrap(~side, scales = "free") +xlab("L2 distance") + ylab("Count")+
  custom_theme

p2 <- ggarrange(g5,g6,g7,g8,nrow = 2,ncol = 2,labels = c("a","b","c","d"))
ggsave(p2,file="../figures/distribution_3.pdf",width=8, height=8,scale = 0.7)


## ---------------------------------------------------------------------------

coef <- read.csv(file = "../data/KeyPlasma/assoc_mat.csv",header = T, row.names = 1, check.names = F)
mat_bin <- as.matrix(coef)
mat_bin[abs(mat_bin)<4] <- 0
mat_bin[mat_bin!=0] <- 1

row_deg <- rowSums(mat_bin)
col_deg <- colSums(mat_bin)

df_deg <- rbind(
  data.frame(degree = row_deg, side = "row"),
  data.frame(degree = col_deg, side = "column")
)

g5 <- ggplot(df_deg[df_deg$side=="column",], aes(x = degree)) +
  geom_histogram(bins = 30, fill="#5495CF",color="white") +
  facet_wrap(~side, scales = "free") + xlab("Degree") + ylab("Count")+
  custom_theme


g6 <- ggplot(df_deg[df_deg$side=="row",], aes(x = degree)) +
  geom_histogram(bins = 30, fill="#9C3F5D",color="white") +
  facet_wrap(~side, scales = "free") +xlab("Degree") + ylab("Count")+
  custom_theme


pdf("../figures/incident_matrix_4.pdf",width = 8.5,height = 4)
Heatmap(
  mat_bin[rowSums(mat_bin)!=0,],
  col = colorRamp2(
    c(min(mat_bin), 0, max(mat_bin)),
    c("white", "white", "#ca0020")
  ),
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  use_raster = FALSE
)
dev.off()

print(nested(mat_bin[rowSums(mat_bin)!=0,], method = "NODF"))
# 38.41459

row_norms <- sqrt(rowSums(coef^2))
row_norms[row_norms == 0] <- 1
coef_row_norm <- coef / row_norms
row_dist <- as.matrix(dist(coef_row_norm, method = "euclidean"))

col_norms <- sqrt(colSums(coef^2))
col_norms[col_norms == 0] <- 1
coef_col_norm <- t(t(coef) / col_norms)
col_dist <- as.matrix(dist(t(coef_col_norm), method = "euclidean"))
row_dist <- reshape2::melt(row_dist)
col_dist <- reshape2::melt(col_dist)


df_l2 <- rbind(
  data.frame(l2 = row_dist$value, side = "row"),
  data.frame(l2 = col_dist$value, side = "column")
)

g7 <- ggplot(df_l2[df_l2$side=="column",], aes(x = l2)) +
  geom_histogram(bins = 30, fill="#5495CF",color="white") +
  facet_wrap(~side, scales = "free") +xlab("L2 distance") + ylab("Count")+
  custom_theme


g8 <- ggplot(df_l2[df_l2$side=="row",], aes(x = l2)) +
  geom_histogram(bins = 30, fill="#9C3F5D",color="white") +
  facet_wrap(~side, scales = "free") +xlab("L2 distance") + ylab("Count")+
  custom_theme

p2 <- ggarrange(g5,g6,g7,g8,nrow = 2,ncol = 2,labels = c("a","b","c","d"))
ggsave(p2,file="../figures/distribution_4.pdf",width=8, height=8,scale = 0.7)
