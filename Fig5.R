############################################################
# gut microbiome-plasma metabilites impact pipeline
# output:
#   1. Figure 6: classfication
#   2. Figure S10: feature combination


library(ggpubr)
library(ggridges)
library(reshape2)
library(ggrepel)
library(stringr)
library(dplyr)
library(tidyr)
library(glmnet)
library(pROC)


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

cardiometabolic_diseases <- c(
  "ACVD",
  "CAD",
  "T2D"
)

# binary classification: CVD vs healthy
disease_label$cvd <- ifelse(disease_label$disease_main %in% cardiometabolic_diseases, 1,
                            ifelse(disease_label$disease_main == "healthy", 0, NA))

disease_label <- disease_label[!is.na(disease_label$cvd), ]

N_train <- read.csv(file = "../data/KeyPlasma/P_train.csv",header = F)
Z_train <- read.csv(file = "../data/KeyPlasma/Z_train_header.csv",header = T)
meta_train <- read.csv(file = "../results/KeyPlasma/meta_train.csv",header = T, row.names = 1, sep = " ")

N_test <- read.csv(file = "../data/KeyPlasma/P_test.csv",header = F)
Z_test <- read.csv(file = "../data/KeyPlasma/Z_test_header.csv",header = T)
meta_test <- read.csv(file = "../results/KeyPlasma/meta_test.csv",header = T, row.names = 1, sep = " ")

convert_meta_to_matrix <- function(meta_df, value_col = "euc") {
  df <- meta_df %>%
    dplyr::select(
      genus = V4,
      sample = V5,
      value = all_of(value_col)
    )
  
  mat <- df %>%
    group_by(sample, genus) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = genus,
      values_from = value,
      values_fill = 0
    ) %>%
    as.data.frame()
  
  rownames(mat) <- mat$sample
  mat$sample <- NULL
  as.matrix(mat)
}

X_micro_healthy <- t(as.matrix(N_train))
X_micro_disease <- t(as.matrix(N_test))

rownames(X_micro_healthy) <- colnames(Z_train)
rownames(X_micro_disease) <- colnames(Z_test)

colnames(X_micro_healthy) <- rownames(Z_train)
colnames(X_micro_disease) <- rownames(Z_test)

X_impact_healthy <- convert_meta_to_matrix(meta_train, "euc")
X_impact_disease <- convert_meta_to_matrix(meta_test, "euc")

X_null_healthy <- convert_meta_to_matrix(meta_train, "euc_null")
X_null_disease <- convert_meta_to_matrix(meta_test, "euc_null")

X_rand_healthy <- convert_meta_to_matrix(meta_train, "euc_rand")
X_rand_disease <- convert_meta_to_matrix(meta_test, "euc_rand")

align_two_matrices <- function(mat1, mat2) {
  all_cols <- union(colnames(mat1), colnames(mat2))
  
  add_missing <- function(mat, all_cols) {
    miss <- setdiff(all_cols, colnames(mat))
    if (length(miss) > 0) {
      zero_block <- matrix(0, nrow = nrow(mat), ncol = length(miss))
      colnames(zero_block) <- miss
      rownames(zero_block) <- rownames(mat)
      mat <- cbind(mat, zero_block)
    }
    mat[, all_cols, drop = FALSE]
  }
  
  mat1 <- add_missing(mat1, all_cols)
  mat2 <- add_missing(mat2, all_cols)
  
  list(mat1 = mat1, mat2 = mat2)
}

tmp <- align_two_matrices(X_micro_healthy, X_micro_disease)
X_micro_healthy <- tmp$mat1
X_micro_disease <- tmp$mat2

tmp <- align_two_matrices(X_impact_healthy, X_impact_disease)
X_impact_healthy <- tmp$mat1
X_impact_disease <- tmp$mat2

tmp <- align_two_matrices(X_null_healthy, X_null_disease)
X_null_healthy <- tmp$mat1
X_null_disease <- tmp$mat2

tmp <- align_two_matrices(X_rand_healthy, X_rand_disease)
X_rand_healthy <- tmp$mat1
X_rand_disease <- tmp$mat2

X_micro <- rbind(X_micro_healthy, X_micro_disease)
X_impact <- rbind(X_impact_healthy, X_impact_disease)
X_null <- rbind(X_null_healthy, X_null_disease)
X_rand <- rbind(X_rand_healthy, X_rand_disease)

X_micro_null <- cbind(X_micro, X_null)
X_micro_impact <- cbind(X_micro, X_impact)
X_delta <- X_impact - X_null
X_micro_delta <- cbind(X_micro, X_delta)


commpn <- (intersect(rownames(disease_label),rownames(X_micro_null)))
X_micro_null <- X_micro_null[commpn,]
X_micro_impact <- X_micro_impact[commpn,]
X_delta <- X_delta[commpn,]
X_micro_delta <- X_micro_delta[commpn,]

print(dim(X_micro_null))
print(dim(X_micro_impact))
print(dim(X_delta))
print(dim(X_micro_delta))

coef <- read.csv("../data/KeyPlasma/assoc_mat.csv", header = TRUE, row.names = 1, check.names = FALSE)
coef <- as.matrix(coef)

X_metab <- X_micro %*% coef

X_micro_metab <- cbind(X_micro, X_metab)

y <- disease_label$cvd[match(rownames(X_micro_null),rownames(disease_label))]
study <- disease_label$study_name[match(rownames(X_micro_null),rownames(disease_label))]


library(ranger)
library(pROC)

run_cv_auc_balanced_rf <- function(X, y, nfolds = 5, nrepeats = 20,
                                   seed = 123,
                                   balance_train = TRUE,
                                   num.trees = 500, mtry = NULL,
                                   min.node.size = 5,
                                   importance = "none") {
  set.seed(seed)
  aucs <- c()
  
  X <- as.data.frame(X)
  y <- as.numeric(y)
  
  if (!all(sort(unique(y)) == c(0, 1))) {
    stop("y must be coded as 0/1")
  }
  
  make_stratified_folds <- function(y, nfolds) {
    fold_id <- integer(length(y))
    
    idx0 <- which(y == 0)
    idx1 <- which(y == 1)
    
    idx0 <- sample(idx0)
    idx1 <- sample(idx1)
    
    fold_id[idx0] <- sample(rep(seq_len(nfolds), length.out = length(idx0)))
    fold_id[idx1] <- sample(rep(seq_len(nfolds), length.out = length(idx1)))
    
    fold_id
  }
  
  for (r in seq_len(nrepeats)) {
    fold_id <- make_stratified_folds(y, nfolds)
    
    for (k in seq_len(nfolds)) {
      test_idx <- which(fold_id == k)
      train_idx <- setdiff(seq_along(y), test_idx)
      
      X_train <- X[train_idx, , drop = FALSE]
      X_test  <- X[test_idx, , drop = FALSE]
      y_train <- y[train_idx]
      y_test  <- y[test_idx]
      
      if (balance_train) {
        idx0 <- which(y_train == 0)
        idx1 <- which(y_train == 1)
        
        n_min <- min(length(idx0), length(idx1))
        
        idx0_sub <- sample(idx0, n_min)
        idx1_sub <- sample(idx1, n_min)
        
        keep_idx <- sample(c(idx0_sub, idx1_sub))
        
        X_train <- X_train[keep_idx, , drop = FALSE]
        y_train <- y_train[keep_idx]
        
        idx0 <- which(y_test == 0)
        idx1 <- which(y_test == 1)
        
        n_min <- min(length(idx0), length(idx1))
        
        idx0_sub <- sample(idx0, n_min)
        idx1_sub <- sample(idx1, n_min)
        
        keep_idx <- sample(c(idx0_sub, idx1_sub))
        
        X_test <- X_test[keep_idx, , drop = FALSE]
        y_test <- y_test[keep_idx]
        
      }
      
      keep <- vapply(X_train, function(col) {
        s <- sd(col, na.rm = TRUE)
        is.finite(s) && s > 0
      }, logical(1))
      
      X_train <- X_train[, keep, drop = FALSE]
      X_test  <- X_test[, keep, drop = FALSE]
      
      dat_train <- data.frame(y = factor(y_train, levels = c(0, 1)), X_train)
      dat_test  <- data.frame(y = factor(y_test, levels = c(0, 1)), X_test)
      
      if (is.null(mtry)) {
        mtry_use <- max(1, floor(sqrt(ncol(X_train))))
      } else {
        mtry_use <- min(mtry, ncol(X_train))
      }
      
      fit <- ranger::ranger(
        y ~ .,
        data = dat_train,
        probability = TRUE,
        num.trees = num.trees,
        mtry = mtry_use,
        min.node.size = min.node.size,
        importance = importance,
        seed = seed + r * 1000 + k
      )
      pred <- predict(fit, data = dat_test)$predictions[, "1"]
      auc_val <- as.numeric(pROC::auc(y_test, pred))
      aucs <- c(aucs, auc_val)
    }
  }
  
  data.frame(mean_auc = aucs)
}

run_loso_auc_rf <- function(X, y, study,
                            balance_train = FALSE,
                            num.trees = 500,
                            seed = 123) {
  
  X <- as.data.frame(X)
  y <- as.numeric(y)
  study <- as.character(study)
  
  aucs <- c()
  study_res <- list()
  effect_study <- c("JieZ_2017","MetaCardis_2020_a","QinJ_2012","SankaranarayananK_2015")
  
  for (s in effect_study) {
    
    test_idx <- which(study == s)
    train_idx <- which(study != s)
    
    X_train <- X[train_idx,]
    X_test  <- X[test_idx,]
    
    y_train <- y[train_idx]
    y_test  <- y[test_idx]
    
    X_train <- as.data.frame(X_train)
    X_test  <- as.data.frame(X_test)
    
    bad_train <- vapply(X_train, function(x) any(!is.finite(as.numeric(x))), logical(1))
    bad_test  <- vapply(X_test,  function(x) any(!is.finite(as.numeric(x))), logical(1))
    
    keep <- !(bad_train | bad_test)
    
    X_train <- X_train[, keep, drop = FALSE]
    X_test  <- X_test[, keep, drop = FALSE]
    
    # optional balancing
    if (balance_train) {
      
      idx0 <- which(y_train == 0)
      idx1 <- which(y_train == 1)
      
      n_min <- min(length(idx0), length(idx1))
      
      keep <- c(
        sample(idx0, n_min),
        sample(idx1, n_min)
      )
      
      X_train <- X_train[keep,]
      y_train <- y_train[keep]
    }
    
    dat_train <- data.frame(
      y = factor(y_train, levels = c(0,1)),
      X_train
    )
    
    dat_test <- data.frame(
      y = factor(y_test, levels = c(0,1)),
      X_test
    )
    
    fit <- ranger::ranger(
      y ~ .,
      data = dat_train,
      probability = TRUE,
      num.trees = 200,
      mtry = max(1, floor(sqrt(ncol(X_train)))),
      min.node.size = 10,
      num.threads = 1,
      seed = seed
    )
    
    pred <- predict(fit, dat_test)$predictions[,1]
    
    auc_val <- as.numeric(
      pROC::auc(y_test, pred)
    )
    
    aucs <- c(aucs, auc_val)
    
    study_res[[s]] <- auc_val
  }
  
  data.frame(
    study = names(study_res),
    auc = unlist(study_res)
  )
}

X_micro_covar <- cbind(X_micro[rownames(X_micro_null), ])
X_metab_covar <- cbind(X_metab[rownames(X_micro_null), ])
X_micro_null_covar <- cbind(X_micro_null)
X_micro_impact_covar <- cbind(X_micro_impact)
X_rand_covar <- cbind(X_rand[rownames(X_micro_null), ])


# res_micro <- run_cv_auc_balanced_rf(X_micro_covar, y)
# res_metab <- run_cv_auc_balanced_rf(X_metab_covar, y)
# res_micro_null <- run_cv_auc_balanced_rf(X_micro_null_covar, y)
# res_micro_impact <- run_cv_auc_balanced_rf(X_micro_impact_covar, y)
# res_rand <- run_cv_auc_balanced_rf(X_rand_covar, y)
# 


loso_micro <- run_loso_auc_rf(
  X_micro_covar,
  y,
  covar_df$study_name
)

loso_metab <- run_loso_auc_rf(
  X_metab_covar,
  y,
  covar_df$study_name
)

loso_impact <- run_loso_auc_rf(
  X_micro_impact_covar,
  y,
  covar_df$study_name
)

loso_null <- run_loso_auc_rf(
  X_micro_null_covar,
  y,
  covar_df$study_name
)

loso_rand <- run_loso_auc_rf(
  X_rand_covar,
  y,
  covar_df$study_name
)


results <- rbind(
  data.frame(model = "Microbiome", loso_micro),
  data.frame(model = "Predicted metabolites", loso_metab),
  data.frame(model = "Microbiome + null impact", loso_null),
  data.frame(model = "Microbiome + impact", loso_impact),
  data.frame(model = "Random impact", loso_rand)
)


g1 <- ggplot(data = results[!results$model=="Microbiome + null impact",], aes(model,auc, fill=model)) +
  coord_flip() + 
  scale_fill_manual(values = c("#80A1C1", "#C94277", "#EEE3AB", "#274C77", "#5E8C61")) + 
  geom_boxplot(linewidth=0.3) + 
  ylab("AUROC") + xlab("") + 
  stat_compare_means(
    comparisons = list(
      c("Microbiome", "Microbiome + impact"),
      c("Microbiome", "Predicted metabolites"),
      c("Random impact", "Microbiome + impact"),
      c("Random impact", "Predicted metabolites")
    ),
    method = "t.test",
    label = "p.signif",
    paired = T
  ) +
  custom_theme + theme(legend.position = "none")


library(glmnet)

run_assoc_lm <- function(X, y, study, min_nonzero = 10) {
  X <- as.data.frame(X)
  y <- as.numeric(y)
  study <- factor(study)
  
  res <- lapply(colnames(X), function(feat) {
    x <- X[[feat]]
    
    # filter: remove near-constant or too sparse features
    if (sd(x, na.rm = TRUE) == 0) {
      return(NULL)
    }
    
    if (sum(x != 0, na.rm = TRUE) < min_nonzero) {
      return(NULL)
    }
    
    df <- data.frame(y = y, x = scale(x), study = study)
    
    fit <- tryCatch(
      glm(y ~ x + study, data = df, family = "binomial"),
      warning = function(w) NULL,
      error = function(e) NULL
    )
    
    if (is.null(fit)) return(NULL)
    
    s <- summary(fit)$coefficients
    
    if (!"x" %in% rownames(s)) return(NULL)
    
    data.frame(
      feature = feat,
      beta = s["x", "Estimate"],
      p = s["x", "Pr(>|z|)"]
    )
  })
  
  res <- do.call(rbind, res)
  res$padj <- p.adjust(res$p, method = "BH")
  res
}

study <- disease_label$study_name[
  match(rownames(X_micro[rownames(X_micro_null),]), rownames(disease_label))
]

res_micro  <- run_assoc_lm(X_micro[rownames(X_micro_null),], y, study)
res_impact <- run_assoc_lm(X_impact[rownames(X_micro_null),], y, study)

res_micro$model <- "Microbiome"
res_impact$model <- "Impact"

res_all <- rbind(
  res_micro[, c("feature", "beta", "p", "padj", "model")],
  res_impact[, c("feature", "beta", "p", "padj", "model")]
)


res_all$logp <- -log10(res_all$padj)
res_all$abs_beta <- abs(res_all$beta)
res_all$logp[is.infinite(res_all$logp)] <- max(res_all$logp[is.finite(res_all$logp)], na.rm = TRUE) + 1


res_all <- res_all %>%
  mutate(
    signif = case_when(
      padj < 0.05 & beta > 0 ~ "FDR<0.05 (positive)",
      padj < 0.05 & beta < 0 ~ "FDR<0.05 (negative)",
      TRUE ~ "NS"
    ),
    logp = -log10(padj)
  )

top_features <- res_all %>%
  filter(padj < 0.001) %>%
  group_by(model) %>%
  arrange(padj) %>%
  slice_head(n = 5) %>%   # top 5 per model
  ungroup()


g_volcano <- ggplot(res_all, aes(x = beta, y = logp)) +
  geom_point(aes(color = signif), alpha = 0.6, size = 1) +
  
  # label top features
  geom_text_repel(
    data = top_features,
    aes(label = feature),
    size = 2,
    max.overlaps = Inf
  ) +
  
  facet_wrap(~model, scales = "free_x") +
  
  scale_color_manual(values = c(
    "FDR<0.05 (positive)" = "#d73027",
    "FDR<0.05 (negative)" = "#4575b4",
    "NS" = "grey80"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed") + 
  custom_theme1 +
  xlab("Effect size (beta)") +
  ylab(expression(-log[10](FDR))) +
  theme(legend.position = "top")


gC <- ggplot(res_all, aes(x = padj, color = model)) +
  stat_ecdf(linewidth = 0.8) +
  scale_color_manual(values = c("Microbiome" = "#4575b4", "Impact" = "#d73027")) +
  geom_vline(xintercept = 0.05, linetype = "dashed") +
  scale_x_log10() + 
  custom_theme1 +
  xlab("Adjusted p-value") +
  ylab("Cumulative fraction of features") +
  theme(legend.position = "none")


p1 <- ggarrange(g1, gC, ncol = 2, nrow = 1, labels = c("a", "c"),widths = c(1.0,0.8))
p2 <- ggarrange(p1, g_volcano, ncol = 1, nrow = 2, labels = c("", "b"),heights = c(0.7,1))
ggsave(p2,file="../figures/classfication.pdf",width=7, height=8,scale = 0.8)


write.csv(results,"../results/auroc.csv")
write.csv(res_all,"../results/lm.csv")

