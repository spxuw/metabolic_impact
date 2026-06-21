############################################################
# gut microbiome-plasma metabilites impact pipeline
# output:
#   1. Figure 7: dietary intervention


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

diet <- read_excel("../data/41467_2023_41042_MOESM9_ESM.xlsx",sheet = 2)
microbiome <- read_excel("../data/41467_2023_41042_MOESM9_ESM.xlsx",sheet = 3)
metabolites <- read_excel("../data/41467_2023_41042_MOESM9_ESM.xlsx",sheet = 7, skip = 12)
Cytokines <- read_excel("../data/41467_2023_41042_MOESM9_ESM.xlsx",sheet = 8)


make_key <- function(df) {
  df %>%
    mutate(
      sample_id = paste(`Participant ID`, `Time Point`, sep = "__")
    )
}

microbiome <- make_key(microbiome)
metabolites <- make_key(metabolites)
Cytokines <- make_key(Cytokines)

# keep only matched samples
common_ids <- intersect(microbiome$sample_id, metabolites$sample_id)

microbiome_sub <- microbiome %>%
  filter(sample_id %in% common_ids)

metabolites_sub <- metabolites %>%
  filter(sample_id %in% common_ids)

# keep only participants with both pre and post in both datasets
paired_ids <- microbiome_sub %>%
  dplyr::count(`Participant ID`, `Time Point`) %>%
  select(`Participant ID`, `Time Point`) %>%
  distinct()

participants_keep <- Reduce(
  intersect,
  list(
    microbiome_sub %>%
      dplyr::count(`Participant ID`, `Time Point`) %>%
      pivot_wider(names_from = `Time Point`, values_from = n, values_fill = 0) %>%
      filter(`Pre-intervention` > 0, `Post-intervention` > 0) %>%
      pull(`Participant ID`),
    metabolites_sub %>%
      dplyr::count(`Participant ID`, `Time Point`) %>%
      pivot_wider(names_from = `Time Point`, values_from = n, values_fill = 0) %>%
      filter(`Pre-intervention` > 0, `Post-intervention` > 0) %>%
      pull(`Participant ID`)
  )
)

microbiome_sub <- microbiome_sub %>%
  filter(`Participant ID` %in% participants_keep)

metabolites_sub <- metabolites_sub %>%
  filter(`Participant ID` %in% participants_keep)


# metadata columns
meta_cols <- c("Diet", "Participant ID", "Time Point", "sample_id")

micro_cols <- setdiff(colnames(microbiome_sub), meta_cols)

extract_genus <- function(x) {
  m <- str_match(x, "\\|g__([^|]+)")
  out <- m[,2]
  out[is.na(out)] <- x[is.na(out)]
  out
}

genus_names <- extract_genus(micro_cols)

# collapse species to genus by row-wise sum
micro_genus <- microbiome_sub %>%
  select(all_of(meta_cols)) %>%
  bind_cols(
    as.data.frame(
      rowsum(
        t(as.matrix(microbiome_sub[, micro_cols, drop = FALSE])),
        group = genus_names,
        reorder = FALSE
      ) |> t()
    )
  )

# optional: replace NA with 0
micro_genus[is.na(micro_genus)] <- 0

coef <- as.matrix(coef)

genus_common <- intersect(colnames(micro_genus)[!(colnames(micro_genus) %in% meta_cols)], rownames(coef))

micro_genus_use <- micro_genus %>%
  select(all_of(meta_cols), all_of(genus_common))

coef_use <- coef[genus_common, , drop = FALSE]

micro_pre <- micro_genus_use %>%
  filter(`Time Point` == "Pre-intervention") %>%
  arrange(`Participant ID`)

micro_post <- micro_genus_use %>%
  filter(`Time Point` == "Post-intervention") %>%
  arrange(`Participant ID`)

stopifnot(all(micro_pre$`Participant ID` == micro_post$`Participant ID`))

X_pre <- as.matrix(micro_pre[, genus_common, drop = FALSE])
X_post <- as.matrix(micro_post[, genus_common, drop = FALSE])

rownames(X_pre) <- micro_pre$`Participant ID`
rownames(X_post) <- micro_post$`Participant ID`

met_meta_cols <- c("Diet", "Participant ID", "Time Point", "sample_id")
met_cols <- setdiff(colnames(metabolites_sub), met_meta_cols)

met_pre <- metabolites_sub %>%
  filter(`Time Point` == "Pre-intervention") %>%
  arrange(`Participant ID`)

met_post <- metabolites_sub %>%
  filter(`Time Point` == "Post-intervention") %>%
  arrange(`Participant ID`)

stopifnot(all(met_pre$`Participant ID` == met_post$`Participant ID`))

Y_pre <- as.matrix(met_pre[, met_cols, drop = FALSE])
Y_post <- as.matrix(met_post[, met_cols, drop = FALSE])

rownames(Y_pre) <- met_pre$`Participant ID`
rownames(Y_post) <- met_post$`Participant ID`

# microbiome change
X_delta <- X_post - X_pre

# metabolite change
Y_delta <- Y_post - Y_pre

# ensure no missing
X_delta[is.na(X_delta)] <- 0
Y_delta[is.na(Y_delta)] <- 0

# standardize
X_delta_z1 <- scale(X_delta[micro_pre$Diet=="MED diet",])
Y_delta_z1 <- scale(Y_delta[met_pre$Diet=="MED diet",])

X_delta_z2 <- scale(X_delta[micro_pre$Diet=="PPT diet",])
Y_delta_z2 <- scale(Y_delta[met_pre$Diet=="PPT diet",])

# association matrix (genus x metabolite)
assoc_obs1 <- t(X_delta_z1) %*% Y_delta_z1 / (nrow(X_delta_z1) - 1)
assoc_obs2 <- t(X_delta_z2) %*% Y_delta_z2 / (nrow(X_delta_z2) - 1)

# attach genus
impact_df <- as.data.frame(diff_mat)
impact_df$genus <- meta$V4

impact_null_df <- as.data.frame(diff_mat_null)
impact_null_df$genus <- meta$V4

# aggregate
impact_genus <- impact_df %>%
  group_by(genus) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

impact_null_genus <- impact_null_df %>%
  group_by(genus) %>%
  summarise(across(everything(), mean, na.rm = TRUE))

impact_mat <- as.matrix(impact_genus[,-1])
rownames(impact_mat) <- impact_genus$genus

impact_null <- as.matrix(impact_null_genus[,-1])
rownames(impact_null) <- impact_null_genus$genus
colnames(impact_null) <- colnames(impact_mat)

common_genus <- intersect(rownames(assoc_obs1),rownames(impact_mat))
common_metabolites <- intersect(colnames(assoc_obs1),colnames(impact_mat))

assoc_obs_common1 <- assoc_obs1[common_genus,common_metabolites]
assoc_obs_common2 <- assoc_obs2[common_genus,common_metabolites]

impact_mat_common <- impact_mat[common_genus,common_metabolites]
impact_mat_null_common <- impact_null[common_genus,common_metabolites]


df <- data.frame(
  assoc1 = as.vector(assoc_obs_common1),
  assoc2 = as.vector(assoc_obs_common2),
  impact = as.vector(impact_mat_common),
  impact_null = as.vector(impact_mat_null_common)
)

# remove NA
df <- df %>% filter(!is.na(assoc1),!is.na(assoc2), !is.na(impact), !is.na(impact_null))

# absolute values for ranking
df <- df %>%
  mutate(
    assoc_abs1 = (assoc1),
    assoc_abs2 = (assoc2),
    impact_abs = (impact),
    impact_null_abs = (impact_null)
  )

rank_obs1 <- order(df$assoc_abs1, decreasing = TRUE)
rank_impact <- order(df$impact_abs, decreasing = TRUE)
rank_null <- order(df$impact_null_abs, decreasing = TRUE)

k_grid <- seq(500,3000,100)

overlap_curve <- lapply(k_grid, function(k) {
  top_obs <- rank_obs1[1:k]
  top_imp <- rank_impact[1:k]
  top_null <- rank_null[1:k]
  
  overlap_imp <- length(intersect(top_obs, top_imp))
  overlap_null <- length(intersect(top_obs, top_null))
  
  data.frame(
    k = k,
    overlap_interaction = overlap_imp,
    overlap_null = overlap_null
  )
})

overlap_curve <- do.call(rbind, overlap_curve)


df_plot <- overlap_curve %>%
  pivot_longer(
    cols = c(overlap_interaction, overlap_null),
    names_to = "model",
    values_to = "overlap"
  )

g1 <- ggplot(df_plot, aes(x = k/nrow(df), y = overlap, color = model)) +
  geom_line(size = 0.5) +
  geom_point(size = 2) +
  custom_theme + 
  xlab("Top-n pairs / total pairs") +
  ylab("Overlap") +
  scale_color_manual(values = c(
    "overlap_interaction" = "#d73027",
    "overlap_null" = "grey50"
  ))


rank_obs1 <- order(df$assoc_abs2, decreasing = TRUE)
rank_impact <- order(df$impact_abs, decreasing = TRUE)
rank_null <- order(df$impact_null_abs, decreasing = TRUE)

k_grid <- seq(500,3000,100)

overlap_curve <- lapply(k_grid, function(k) {
  top_obs <- rank_obs1[1:k]
  top_imp <- rank_impact[1:k]
  top_null <- rank_null[1:k]
  
  overlap_imp <- length(intersect(top_obs, top_imp))
  overlap_null <- length(intersect(top_obs, top_null))
  
  data.frame(
    k = k,
    overlap_interaction = overlap_imp,
    overlap_null = overlap_null
  )
})

overlap_curve <- do.call(rbind, overlap_curve)


df_plot <- overlap_curve %>%
  pivot_longer(
    cols = c(overlap_interaction, overlap_null),
    names_to = "model",
    values_to = "overlap"
  )

g11 <- ggplot(df_plot, aes(x = k/nrow(df), y = overlap, color = model)) +
  geom_line(size = 0.5) +
  geom_point(size = 2) +
  custom_theme + 
  xlab("Top-n pairs / total pairs") +
  ylab("Overlap") +
  scale_color_manual(values = c(
    "overlap_interaction" = "#d73027",
    "overlap_null" = "grey50"
  ))

p1 <- ggarrange(g1,g11,nrow = 1,ncol = 2,labels = c("a","b"))
ggsave(p1,file="../figures/overlap_comparison.pdf",width=11, height=4,scale = 0.8)


# ==============================
# Function to extract interaction-specific overlapping pairs
# ==============================
get_interaction_specific_overlap <- function(df_full,
                                             obs_col = "assoc_abs1",
                                             impact_col = "impact_abs",
                                             null_col = "impact_null_abs",
                                             top_n = 3000,
                                             heatmap_top_n = 80) {
  
  df_ranked <- df_full %>%
    mutate(
      pair_id = paste(genus, metabolite, sep = "___"),
      rank_obs = rank(-.data[[obs_col]], ties.method = "first"),
      rank_imp = rank(-.data[[impact_col]], ties.method = "first"),
      rank_null = rank(-.data[[null_col]], ties.method = "first")
    )
  
  top_obs  <- df_ranked %>% filter(rank_obs <= top_n)  %>% pull(pair_id)
  top_imp  <- df_ranked %>% filter(rank_imp <= top_n)  %>% pull(pair_id)
  top_null <- df_ranked %>% filter(rank_null <= top_n) %>% pull(pair_id)
  
  # observed + interaction-aware overlap
  overlap_imp <- intersect(top_obs, top_imp)
  
  # observed + null overlap
  overlap_null <- intersect(top_obs, top_null)
  
  # interaction-aware-specific pairs
  specific_pairs <- setdiff(overlap_imp, overlap_null)
  
  df_specific <- df_ranked %>%
    filter(pair_id %in% specific_pairs) %>%
    mutate(
      combined_rank = rank_obs + rank_imp
    ) %>%
    arrange(combined_rank)
  
  # restrict to manageable number for heatmap
  df_specific_top <- df_specific %>%
    slice_head(n = heatmap_top_n)
  
  list(
    df_specific = df_specific,
    df_specific_top = df_specific_top,
    n_overlap_imp = length(overlap_imp),
    n_overlap_null = length(overlap_null),
    n_specific = length(specific_pairs)
  )
}

df_full <- data.frame(
  genus = rep(rownames(assoc_obs_common1), times = ncol(assoc_obs_common1)),
  metabolite = rep(colnames(assoc_obs_common1), each = nrow(assoc_obs_common1)),
  assoc1 = as.vector(assoc_obs_common1),
  assoc2 = as.vector(assoc_obs_common2),
  impact = as.vector(impact_mat_common),
  impact_null = as.vector(impact_mat_null_common)
) %>%
  mutate(
    assoc_abs1 = assoc1,
    assoc_abs2 = assoc2,
    impact_abs = impact,
    impact_null_abs = impact_null
  ) %>%
  filter(
    !is.na(assoc_abs1),
    !is.na(assoc_abs2),
    !is.na(impact_abs),
    !is.na(impact_null_abs)
  )


res_med <- get_interaction_specific_overlap(
  df_full,
  obs_col = "assoc_abs1",
  impact_col = "impact_abs",
  null_col = "impact_null_abs",
  top_n = 3000,
  heatmap_top_n = 80
)

cat("MED interaction overlap:", res_med$n_overlap_imp, "\n")
cat("MED null overlap:", res_med$n_overlap_null, "\n")
cat("MED interaction-specific:", res_med$n_specific, "\n")

top_overlap_med_specific <- res_med$df_specific_top

mat_df <- top_overlap_med_specific %>%
  select(genus, metabolite, impact) %>%
  pivot_wider(
    names_from = metabolite,
    values_from = impact,
    values_fill = 0
  )

mat <- as.matrix(mat_df[, -1])
rownames(mat) <- mat_df$genus

pdf("../figures/heatmap_interaction_specific_MED_top3000.pdf",
    width = 7, height = 5)
pheatmap(
  mat,
  fontsize = 6,
  main = "MED: observed + interaction-aware, not null"
)
dev.off()


res_ppt <- get_interaction_specific_overlap(
  df_full,
  obs_col = "assoc_abs2",
  impact_col = "impact_abs",
  null_col = "impact_null_abs",
  top_n = 3000,
  heatmap_top_n = 80
)

cat("PPT interaction overlap:", res_ppt$n_overlap_imp, "\n")
cat("PPT null overlap:", res_ppt$n_overlap_null, "\n")
cat("PPT interaction-specific:", res_ppt$n_specific, "\n")

top_overlap_ppt_specific <- res_ppt$df_specific_top


mat_df <- top_overlap_ppt_specific %>%
  select(genus, metabolite, impact) %>%
  pivot_wider(
    names_from = metabolite,
    values_from = impact,
    values_fill = 0
  )

mat <- as.matrix(mat_df[, -1])
rownames(mat) <- mat_df$genus

pdf("../figures/heatmap_interaction_specific_PPT_top3000.pdf",
    width = 7, height = 3.5)
pheatmap(
  mat,
  fontsize = 6,
  main = "PPT: observed + interaction-aware, not null"
)
dev.off()


# assocaition 
N_train <- read.csv(file = "../data/KeyPlasma/P_train.csv",header = F)
Z_train <- read.csv(file = "../data/KeyPlasma/Z_train_header.csv",header = T)
rownames(N_train) <- rownames(Z_train)

# compute variance per genus (row)
row_var <- apply(N_train, 1, var)
N_train_filt <- N_train[row_var > 0, ]

library(Hmisc)

# rcorr expects samples in rows → transpose again
res <- rcorr(t(N_train_filt), type = "spearman")

cor_genus <- res$r
p_genus   <- res$P

target <- unique(top_overlap_med_specific$genus)

cor_genus_sub <- cor_genus[target, target]
p_genus_sub   <- p_genus[target, target]

diag(cor_genus_sub) <- NA
diag(p_genus_sub)   <- NA

p_adj <- matrix(
  p.adjust(p_genus_sub, method = "BH"),
  nrow = nrow(p_genus_sub),
  dimnames = dimnames(p_genus_sub)
)

sig_labels <- matrix("", nrow = nrow(p_adj), ncol = ncol(p_adj))

sig_labels[p_adj < 0.05]  <- "*"
sig_labels[p_adj < 0.01]  <- "**"
sig_labels[p_adj < 0.001] <- "***"

pdf("../figures/heatmap_med_focus.pdf",width = 8,height = 7)
pheatmap(
  cor_genus_sub,
  display_numbers = sig_labels,
  number_color = "black",
  main = "Spearman correlation (FDR-adjusted significance)"
)
dev.off()



library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE76925", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }


library(illuminaHumanv4.db)
library(AnnotationDbi)


probe_ids <- rownames(ex)

anno <- AnnotationDbi::select(
  illuminaHumanv4.db,
  keys = probe_ids,
  columns = c("SYMBOL", "ENTREZID", "GENENAME"),
  keytype = "PROBEID"
)

ex_df <- as.data.frame(ex)
ex_df$PROBEID <- rownames(ex_df)

ex_anno <- merge(
  anno,
  ex_df,
  by.x = "PROBEID",
  by.y = "PROBEID"
)

head(ex_anno[, 1:5])

library(dplyr)

ex_gene <- ex_anno %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  group_by(SYMBOL) %>%
  summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

ex_gene_mat <- as.matrix(ex_gene[, -1])
rownames(ex_gene_mat) <- ex_gene$SYMBOL

write.table(ex,file = "GSE76925_expr.csv", row.names = T, col.names = T, sep = ",")

ex_anno_valid <- ex_anno[complete.cases(ex_anno[,1:2]),]
write.table(ex_anno_valid[,1:2],file = "ex_anno.csv", row.names = F, col.names = T, sep = ",")

write.table(ex_gene_mat,file = "GSE76925_expr_genelevel.csv", row.names = T, col.names = T, sep = ",")


meta <- gset@phenoData@data
meta <- meta[colnames(ex_gene_mat),]

write.table(meta,file = "GSE76925_meta.csv", row.names = T, col.names = T, sep = ",")
