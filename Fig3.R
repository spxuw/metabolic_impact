############################################################
# gut microbiome-plasma metabilites impact pipeline
# output:
#   1. Figure 3: top taxa impact distribution
#   2. Figure 4: differential analysis


library(ggpubr)
library(ggridges)
library(reshape2)
library(ggrepel)


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

meta_train <- read.csv(file = "../results/KeyPlasma/meta_train.csv",sep = " ")

Prevlance = as.data.frame(table(meta_train$V4))
Prevlance = Prevlance[order(Prevlance$Freq,decreasing = T),]

Prevlance = Prevlance[Prevlance$Freq>length(unique(meta_train$V3))*0.05,]
meta_train = meta_train[meta_train$V4%in%Prevlance$Var1,]

keystoness_median_train = aggregate(meta_train$euc, list(meta_train$V4), FUN=median)
keystoness_median_train = keystoness_median_train[order(keystoness_median_train$x),]


g1 = ggplot(data=meta_train[meta_train$V4%in%keystoness_median_train$Group.1[101:110],],aes(x=euc,y=V4))+
  geom_density_ridges_gradient(fill="#D77186")+
  scale_x_continuous(expand = c(0, 0),limits = c(0,2))+
  xlab('')+ylab('')+custom_theme


g2 = ggplot(data=meta_train[meta_train$V4%in%keystoness_median_train$Group.1[1:10],],aes(x=euc,y=V4))+
  geom_density_ridges_gradient(fill="#3388bd")+
  scale_x_continuous(expand = c(0, 0),limits = c(0,2))+
  xlab('')+ylab('')+custom_theme


meta_test <- read.csv(file = "../results/KeyPlasma/meta_test.csv",sep = " ")

Prevlance = as.data.frame(table(meta_test$V4))
Prevlance = Prevlance[order(Prevlance$Freq,decreasing = T),]

Prevlance = Prevlance[Prevlance$Freq>length(unique(meta_test$V3))*0.05,]
meta_test = meta_test[meta_test$V4%in%Prevlance$Var1,]

keystoness_median_test = aggregate(meta_test$euc, list(meta_test$V4), FUN=median)
keystoness_median_test = keystoness_median_test[order(keystoness_median_test$x),]


g3 = ggplot(data=meta_test[meta_test$V4%in%keystoness_median_test$Group.1[102:111],],aes(x=euc,y=V4))+
  geom_density_ridges_gradient(fill="#D77186")+
  scale_x_continuous(expand = c(0, 0),limits = c(0,2))+
  xlab('')+ylab('')+custom_theme


g4 = ggplot(data=meta_test[meta_test$V4%in%keystoness_median_test$Group.1[1:10],],aes(x=euc,y=V4))+
  geom_density_ridges_gradient(fill="#3388bd")+
  scale_x_continuous(expand = c(0, 0),limits = c(0,2))+
  xlab('')+ylab('')+custom_theme


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

diff_mat <- t(vapply(seq_len(nrow(meta)), function(i) {
  m_col <- N1[meta$V3[i], ]
  n_row <- M1[i, ]
  n_row - m_col
}, numeric(ncol(M1))))

M <- read.csv(file = "../results/KeyPlasma/qtst2.csv",header = F)
N <- read.csv(file = "../data/KeyPlasma/P_test.csv",header = F)
meta_test <- read.csv(file = "../data/KeyPlasma/meta_testr.csv",header = F)

M1 <- as.matrix(M) %*% as.matrix(coef)
N1 <- as.matrix(t(N)) %*% as.matrix(coef)

M1 <- scale(M1, center = mu, scale = sdv)
N1 <- scale(N1, center = mu, scale = sdv)

diff_mat_test <- t(vapply(seq_len(nrow(meta_test)), function(i) {
  m_col <- N1[meta_test$V3[i], ]
  n_row <- M1[i, ]
  n_row - m_col
}, numeric(ncol(M1))))


annotation_met <- read.csv("../data/Supplementary Table 2.tsv",
                           sep = "\t",
                           check.names = FALSE)

annot <- annotation_met %>%
  transmute(
    metabolite = Metabolite,
    class = `Metabolite class`,
    subclass = `Metabolite subclass`
  ) %>%
  distinct()

# Make sure columns of diff_mat are metabolites
colnames(diff_mat) <- colnames(M1)
colnames(diff_mat_test) <- colnames(M1)

common_met <- intersect(colnames(diff_mat), annot$metabolite)

annot_use <- annot %>%
  arrange(match(metabolite, colnames(M1)))

# squared contribution per metabolite
sq_diff <- diff_mat^2
sq_diff_test <- diff_mat_test^2


get_top_contributing_metabolites <- function(diff_mat, meta_df, taxa_use,
                                             top_n = 100) {
  
  idx <- which(meta_df$V4 %in% taxa_use)
  
  contrib <- colMeans(diff_mat[idx, , drop = FALSE]^2, na.rm = TRUE)
  
  data.frame(
    metabolite = names(contrib),
    contribution = as.numeric(contrib)
  ) %>%
    arrange(desc(contribution)) %>%
    slice_head(n = top_n)
}


# ------------------------------------------------------------
# 2. Enrichment function
# ------------------------------------------------------------
run_metabolite_enrichment <- function(top_metabs, background_metabs,
                                      annot_use,
                                      level = c("subclass", "class")) {
  
  level <- match.arg(level)
  
  annot_bg <- annot_use %>%
    filter(metabolite %in% background_metabs) %>%
    filter(!is.na(.data[[level]]), .data[[level]] != "") %>%
    distinct(metabolite, .data[[level]])
  
  annot_fg <- annot_use %>%
    filter(metabolite %in% top_metabs) %>%
    filter(!is.na(.data[[level]]), .data[[level]] != "") %>%
    distinct(metabolite, .data[[level]])
  
  terms <- unique(annot_bg[[level]])
  
  res <- lapply(terms, function(term) {
    
    fg_in <- sum(annot_fg[[level]] == term, na.rm = TRUE)
    fg_out <- nrow(annot_fg) - fg_in
    
    bg_in <- sum(annot_bg[[level]] == term, na.rm = TRUE) - fg_in
    bg_out <- nrow(annot_bg) - nrow(annot_fg) - bg_in
    
    if (fg_in == 0) return(NULL)
    
    ft <- fisher.test(matrix(c(fg_in, fg_out, bg_in, bg_out), nrow = 2))
    
    data.frame(
      term = term,
      overlap = fg_in,
      foreground = nrow(annot_fg),
      background = nrow(annot_bg),
      odds_ratio = unname(ft$estimate),
      p = ft$p.value
    )
  })
  
  bind_rows(res) %>%
    mutate(padj = p.adjust(p, method = "BH")) %>%
    arrange(p)
}


# ------------------------------------------------------------
# 3. Define top/bottom taxa
# ------------------------------------------------------------
top_taxa_train <- keystoness_median_train$Group.1[101:110]
bottom_taxa_train <- keystoness_median_train$Group.1[1:10]

top_taxa_test <- keystoness_median_test$Group.1[102:111]
bottom_taxa_test <- keystoness_median_test$Group.1[1:10]


# ------------------------------------------------------------
# 4. Get top-contributing metabolites
# ------------------------------------------------------------
background_train <- intersect(colnames(diff_mat), annot_use$metabolite)
background_test  <- intersect(colnames(diff_mat_test), annot_use$metabolite)

top_met_train <- get_top_contributing_metabolites(
  diff_mat, meta, top_taxa_train, top_n = 100
)

bottom_met_train <- get_top_contributing_metabolites(
  diff_mat, meta, bottom_taxa_train, top_n = 100
)

top_met_test <- get_top_contributing_metabolites(
  diff_mat_test, meta_test, top_taxa_test, top_n = 100
)

bottom_met_test <- get_top_contributing_metabolites(
  diff_mat_test, meta_test, bottom_taxa_test, top_n = 100
)


# ------------------------------------------------------------
# 5. Run subclass enrichment
# ------------------------------------------------------------
enrich_train_top <- run_metabolite_enrichment(
  top_metabs = top_met_train$metabolite,
  background_metabs = background_train,
  annot_use = annot_use,
  level = "subclass"
) %>% mutate(group = "Healthy top-impact")

enrich_train_bottom <- run_metabolite_enrichment(
  top_metabs = bottom_met_train$metabolite,
  background_metabs = background_train,
  annot_use = annot_use,
  level = "subclass"
) %>% mutate(group = "Healthy bottom-impact")

enrich_test_top <- run_metabolite_enrichment(
  top_metabs = top_met_test$metabolite,
  background_metabs = background_test,
  annot_use = annot_use,
  level = "subclass"
) %>% mutate(group = "Disease top-impact")

enrich_test_bottom <- run_metabolite_enrichment(
  top_metabs = bottom_met_test$metabolite,
  background_metabs = background_test,
  annot_use = annot_use,
  level = "subclass"
) %>% mutate(group = "Disease bottom-impact")

enrich_train_top <- enrich_train_top %>%
  filter(overlap >= 2)

enrich_test_top <- enrich_test_top %>%
  filter(overlap >= 2)


enrich_all <- bind_rows(
  enrich_train_top,
  enrich_test_top,
  enrich_train_bottom,
  enrich_test_bottom
)


plot_enrichment <- function(enrich_df, n_terms = 10, title = "") {
  
  df_plot <- enrich_df %>%
    arrange(p) %>%
    slice_head(n = n_terms)
  
  ggplot(df_plot,
         aes(x = odds_ratio,
             y = reorder(term, odds_ratio))) +
    geom_point(aes(size = overlap, color = -log10(p))) +
    scale_color_viridis_c(option = "D") +
    custom_theme1 +
    xlab("Odds ratio") +
    ylab("") +
    labs(
      title = title,
      size = "Overlap",
      color = "-log10(p)"
    )
}

g_enrich_train_top <- plot_enrichment(enrich_train_top, title = "Healthy top-impact")
g_enrich_test_top <- plot_enrichment(enrich_test_top, title = "Disease top-impact")
g_enrich_train_bottom <- plot_enrichment(enrich_train_bottom, title = "Healthy bottom-impact")
g_enrich_test_bottom <- plot_enrichment(enrich_test_bottom, title = "Disease bottom-impact")


p_enrich <- ggarrange(
  g_enrich_train_top,
  g_enrich_test_top,
  nrow = 2,
  ncol = 1,
  labels = c("e", "f")
)


p1 <- ggarrange(g1,g3,g2,g4,nrow = 2,ncol = 2,labels = c("a","b","c","d"))
p2 <- ggarrange(p1,p_enrich,nrow = 1,ncol = 2,align = "hv")
ggsave(p2,file="../figures/top_taxa.pdf",width=13, height=7,scale = 0.8)



meta_train$dataset <- "train"
meta_test$dataset  <- "test"

meta_all <- rbind(meta_train, meta_test)

df_stat <- meta_all %>%
  group_by(V4) %>%
  summarise(
    n_train = sum(dataset == "train"),
    n_test  = sum(dataset == "test"),
    
    med_train = median(euc[dataset == "train"], na.rm = TRUE),
    med_test  = median(euc[dataset == "test"], na.rm = TRUE),
    
    diff = med_test - med_train,
    
    pval = tryCatch(
      wilcox.test(euc ~ dataset)$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    padj = p.adjust(pval, method = "BH"),
    direction = case_when(
      diff > 0 ~ "increase",
      diff < 0 ~ "decrease",
      TRUE ~ "no_change"
    )
  )

top_inc <- df_stat %>%
  filter(direction == "increase") %>%
  arrange(padj, desc(diff)) %>%
  head(10)

top_dec <- df_stat %>%
  filter(direction == "decrease") %>%
  arrange(padj, diff) %>%
  head(10)

df_label <- df_stat %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(10)


g5 <- ggplot(df_stat, aes(x = diff, y = -log10(padj))) +
  geom_point(aes(color = direction), alpha = 0.7) +
  scale_color_manual(values = c("#4575b4","#d73027","white")) + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(
    data = df_label,
    aes(label = V4),
    size = 3,
    max.overlaps = Inf
  ) +
  custom_theme1+ theme(legend.position = "bottom") +
  xlab("Median difference (disease - healthy)") +
  ylab("-log10 adjusted p-value")


# read microbiome abundance matrices if needed
N_train <- read.csv("../data/KeyPlasma/P_train.csv", header = FALSE, check.names = FALSE)
N_test  <- read.csv("../data/KeyPlasma/P_test.csv",  header = FALSE, check.names = FALSE)
N_train_header <- read.csv("../data/KeyPlasma/Z_train_header.csv", header = T, check.names = FALSE)


N_train <- as.matrix(N_train)
N_test  <- as.matrix(N_test)

# simple abundance effect size:
# log10 median abundance difference with pseudocount
pseudo <- 1e-6

abund_stat <- data.frame(
  V4 = rownames(N_train_header),
  med_abund_train = apply(N_train, 1, mean, na.rm = TRUE),
  med_abund_test  = apply(N_test,  1, mean, na.rm = TRUE),
  mean_abund_train = rowMeans(N_train, na.rm = TRUE),
  mean_abund_test  = rowMeans(N_test,  na.rm = TRUE)
) %>%
  mutate(
    log10_med_train = log10(med_abund_train + pseudo),
    log10_med_test  = log10(med_abund_test + pseudo),
    abund_diff = log10_med_test - log10_med_train
  )

df_plot <- df_stat %>%
  dplyr::inner_join(abund_stat, by = "V4") %>%
  mutate(
    sig = case_when(
      !is.na(padj) & padj < 0.05 & diff > 0 ~ "Impact higher in disease",
      !is.na(padj) & padj < 0.05 & diff < 0 ~ "Impact lower in disease",
      TRUE ~ "NS"
    )
  )

df_label2 <- df_plot %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(10)


g6_new <- ggplot(df_plot, aes(x = diff, y = abund_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(aes(color = sig), alpha = 0.75, size = 2) +
  stat_cor(method = "spearman")+
  geom_text_repel(
    data = df_label2,
    aes(label = V4),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c(
    "Impact higher in disease" = "#d73027",
    "Impact lower in disease" = "#4575b4",
    "NS" = "grey75"
  )) +
  custom_theme1 +
  theme(legend.position = "bottom") +
  xlab("Median impact difference (disease - healthy)") +
  ylab("Log10 mean abundance difference (disease - healthy)")

df_plot_sub <- df_plot[complete.cases(df_plot),]
print(cor(df_plot_sub$diff,df_plot_sub$abund_diff, method = "spearman")) # 0.2564814

p2 <- ggarrange(g5,g6_new,nrow = 1,ncol = 2,labels = c("a","b"))
ggsave(p2,file="../figures/diff_taxa.pdf",width=9, height=5.3,scale = 0.8)




