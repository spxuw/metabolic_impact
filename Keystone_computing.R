############################################################
# gut microbiome-plasma metabilites impact pipeline
# output:
#   1. Compute impact
#   2. Figure S5: Null versus real impact
#   3. Figure S6: rand versus real impact
#   4. Figure S7: L2 versus cosine correlation


library(vegan)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lsa)

setwd("/Users/xuwenwang/Library/CloudStorage/Dropbox/Projects/KeyPlasma/code")

M <- read.csv(file = "../results/KeyPlasma/qtst1.csv",header = F)
N <- read.csv(file = "../data/KeyPlasma/P_train.csv",header = F)
meta <- read.csv(file = "../data/KeyPlasma/meta_trainr.csv",header = F)
coef <- read.csv(file = "../data/KeyPlasma/assoc_mat.csv",header = T, row.names = 1, check.names = F)

M1 <- as.matrix(M) %*% as.matrix(coef)
N1 <- as.matrix(t(N)) %*% as.matrix(coef)

mu <- colMeans(N1, na.rm = TRUE)
sdv <- apply(N1, 2, sd, na.rm = TRUE)
sdv[sdv == 0] <- 1

M1_scaled <- scale(M1, center = mu, scale = sdv)
N1_scaled <- scale(N1, center = mu, scale = sdv)

print(dim(meta))
print(dim(M1_scaled))
print(dim(N1_scaled))

eu_dist <- vapply(seq_len(nrow(meta)), function(i) {
  m_col <- N1_scaled[meta$V3[i],]
  n_row <- M1_scaled[i, ]
  
  m_col <- m_col / sqrt(sum(m_col^2))
  n_row <- n_row / sqrt(sum(n_row^2))
  
  vegdist(rbind(n_row, m_col), method = "euclidean")[1]
}, numeric(1))
meta$euc <- eu_dist

cos_dist <- vapply(seq_len(nrow(meta)), function(i) {
  m_col <- N1_scaled[meta$V3[i], ]
  n_row <- M1_scaled[i, ]
  
  1 - cosine(m_col, n_row)
}, numeric(1))
meta$cos <- cos_dist

eu_dist_null <- c()
for (i in 1:nrow(meta))
{
  print(i)
  m_col <- N[,meta$V3[i]]
  m_col[meta$V2[i]] <- 0
  n_row <- N[,meta$V3[i]]
  m_col1 <- as.matrix(t(m_col)) %*% as.matrix(coef)
  n_row1 <- as.matrix(t(n_row)) %*% as.matrix(coef)
  
  m_col1 <- (m_col1 - mu) / sdv
  n_row1 <- (n_row1 - mu) / sdv
  
  m_col1 <- m_col1 / sqrt(sum(m_col1^2))
  n_row1 <- n_row1 / sqrt(sum(n_row1^2))
  
  eu_dist_null = c(eu_dist_null,vegdist(rbind(n_row1, m_col1), method = "euclidean")[1])
}
meta$euc_null <- eu_dist_null

set.seed(123)
coef_rand <- matrix(sample(as.matrix(coef)), nrow = nrow(coef))
M1_rand <- as.matrix(M) %*% as.matrix(coef_rand)
N1_rand <- as.matrix(t(N)) %*% as.matrix(coef_rand)

M1_rand <- scale(M1_rand, center = mu, scale = sdv)
N1_rand <- scale(N1_rand, center = mu, scale = sdv)

eu_dist_rand <- vapply(seq_len(nrow(meta)), function(i) {
  m_col <- N1_rand[meta$V3[i],]
  n_row <- M1_rand[i, ]

  m_col <- m_col / sqrt(sum(m_col^2))
  n_row <- n_row / sqrt(sum(n_row^2))
  
  vegdist(rbind(n_row, m_col), method = "euclidean")[1]
}, numeric(1))
meta$euc_rand <- eu_dist_rand

write.table(meta,file = "../results/KeyPlasma/meta_train.csv")


diff_mat <- t(vapply(seq_len(nrow(meta)), function(i) {
  m_col <- N1_scaled[meta$V3[i], ]
  n_row <- M1_scaled[i, ]
  
  m_col <- (m_col - mu) / sdv
  n_row <- (n_row - mu) / sdv
  
  n_row - m_col
  
}, numeric(ncol(M1))))


vuln <- colMeans(abs(diff_mat), na.rm = TRUE)
vuln_var <- apply(diff_mat, 2, var, na.rm = TRUE)
vuln_q <- apply(abs(diff_mat), 2, quantile, probs = 0.9, na.rm = TRUE)

df_vuln <- data.frame(
  metabolite = colnames(M1),
  mean_abs = vuln,
  var = vuln_var,
  q90 = vuln_q
)

df_vuln <- df_vuln[order(df_vuln$mean_abs, decreasing = TRUE), ]
write.table(df_vuln,file = "../results/KeyPlasma/vuln_train.csv")



M <- read.csv(file = "../results/KeyPlasma/qtst2.csv",header = F)
N <- read.csv(file = "../data/KeyPlasma/P_test.csv",header = F)
meta_test <- read.csv(file = "../data/KeyPlasma/meta_testr.csv",header = F)

M1 <- as.matrix(M) %*% as.matrix(coef)
N1 <- as.matrix(t(N)) %*% as.matrix(coef)

M1_scaled <- scale(M1, center = mu, scale = sdv)
N1_scaled <- scale(N1, center = mu, scale = sdv)

print(dim(meta_test))
print(dim(M1_scaled))
print(dim(N1_scaled))

eu_dist <- vapply(seq_len(nrow(meta_test)), function(i) {
  m_col <- N1_scaled[meta_test$V3[i],]
  n_row <- M1_scaled[i, ]
  
  m_col <- m_col / sqrt(sum(m_col^2))
  n_row <- n_row / sqrt(sum(n_row^2))
  
  vegdist(rbind(n_row, m_col), method = "euclidean")[1]
}, numeric(1))
meta_test$euc <- eu_dist

cos_dist <- vapply(seq_len(nrow(meta_test)), function(i) {
  m_col <- N1_scaled[meta_test$V3[i], ]
  n_row <- M1_scaled[i, ]
  
  1 - cosine(m_col, n_row)
}, numeric(1))
meta_test$cos <- cos_dist


eu_dist_null <- c()
for (i in 1:nrow(meta_test))
{
  print(i)
  m_col <- N[,meta_test$V3[i]]
  m_col[meta_test$V2[i]] <- 0
  n_row <- N[,meta_test$V3[i]]
  m_col1 <- as.matrix(t(m_col)) %*% as.matrix(coef)
  n_row1 <- as.matrix(t(n_row)) %*% as.matrix(coef)
  
  m_col1 <- (m_col1 - mu) / sdv
  n_row1 <- (n_row1 - mu) / sdv
  
  m_col1 <- m_col1 / sqrt(sum(m_col1^2))
  n_row1 <- n_row1 / sqrt(sum(n_row1^2))
  
  eu_dist_null = c(eu_dist_null,vegdist(rbind(n_row1, m_col1), method = "euclidean")[1])
}
meta_test$euc_null <- eu_dist_null

set.seed(123)
coef_rand <- matrix(sample(as.matrix(coef)), nrow = nrow(coef))

M1_rand <- as.matrix(M) %*% as.matrix(coef_rand)
N1_rand <- as.matrix(t(N)) %*% as.matrix(coef_rand)

M1_rand <- scale(M1_rand, center = mu, scale = sdv)
N1_rand <- scale(N1_rand, center = mu, scale = sdv)

eu_dist_rand <- vapply(seq_len(nrow(meta_test)), function(i) {
  m_col <- N1_rand[meta_test$V3[i],]
  n_row <- M1_rand[i, ]
  
  m_col <- m_col / sqrt(sum(m_col^2))
  n_row <- n_row / sqrt(sum(n_row^2))
  
  vegdist(rbind(n_row, m_col), method = "euclidean")[1]
}, numeric(1))
meta_test$euc_rand <- eu_dist_rand


write.table(meta_test,file = "../results/KeyPlasma/meta_test.csv")


diff_mat_test <- t(vapply(seq_len(nrow(meta_test)), function(i) {
  m_col <- N1_scaled[meta_test$V3[i], ]
  n_row <- M1_scaled[i, ]
  
  m_col <- (m_col - mu) / sdv
  n_row <- (n_row - mu) / sdv
  
  n_row - m_col
}, numeric(ncol(M1))))


vuln_test <- colMeans(abs(diff_mat_test), na.rm = TRUE)
vuln_var_test <- apply(diff_mat_test, 2, var, na.rm = TRUE)
vuln_q_test <- apply(abs(diff_mat_test), 2, quantile, probs = 0.9, na.rm = TRUE)

df_vuln_test <- data.frame(
  metabolite = colnames(M1),
  mean_abs = vuln_test,
  var = vuln_var_test,
  q90 = vuln_q_test
)

df_vuln_test <- df_vuln_test[order(df_vuln_test$mean_abs, decreasing = TRUE), ]
write.table(df_vuln_test,file = "../results/KeyPlasma/vuln_test.csv")



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


g1 <- ggplot(data = meta,aes(euc,euc_null)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10") +
  xlim(0,2) + ylim(0,2) + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  stat_cor() + custom_theme + 
  xlab("L2") + ylab("L2 null")
  

g2 <- ggplot(data = meta_test,aes(euc,euc_null)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10") +
  xlim(0,2) + ylim(0,2) + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  stat_cor() + custom_theme + 
  xlab("L2") + ylab("L2 null")


p1 <- ggarrange(g1,g2,nrow = 1,ncol = 2,labels = c("a","b"))
ggsave(p1,file="../figures/l2_cor.pdf",width=9, height=4,scale = 0.8)


dat1 <- data.frame(l2=c(meta$euc,meta$euc_rand),model=c(rep("True",nrow(meta)),rep("Rand",nrow(meta))))

gr1 <- ggplot(data = dat1,aes(model,l2,fill=model)) + 
  geom_boxplot() + stat_compare_means() + 
  scale_fill_manual(values = c("#9467BDFF","#8C564BFF"))+
  custom_theme + xlab("") + ylab("L2 distance")

dat1 <- data.frame(l2=c(meta_test$euc,meta_test$euc_rand),model=c(rep("True",nrow(meta_test)),rep("Rand",nrow(meta_test))))

gr2 <- ggplot(data = dat1,aes(model,l2,fill=model)) + 
  geom_boxplot() +stat_compare_means() + 
  scale_fill_manual(values = c("#9467BDFF","#8C564BFF"))+
  custom_theme + xlab("") + ylab("L2 distance")

pr <- ggarrange(gr1,gr2,nrow = 1,ncol = 2,labels = c("a","b"))
ggsave(pr,file="../figures/rand.pdf",width=8, height=5,scale = 0.8)
  

dat1 <- data.frame(l2=c(meta$euc,meta$euc_null),model=c(rep("True",nrow(meta)),rep("Null",nrow(meta))))

gr1 <- ggplot(data = dat1,aes(model,l2,fill=model)) + 
  geom_boxplot() + stat_compare_means() + 
  scale_fill_manual(values = c("#9467BDFF","#8C564BFF"))+
  custom_theme + xlab("") + ylab("L2 distance")

dat1 <- data.frame(l2=c(meta_test$euc,meta_test$euc_null),model=c(rep("True",nrow(meta_test)),rep("Null",nrow(meta_test))))

gr2 <- ggplot(data = dat1,aes(model,l2,fill=model)) + 
  geom_boxplot() +stat_compare_means() + 
  scale_fill_manual(values = c("#9467BDFF","#8C564BFF"))+
  custom_theme + xlab("") + ylab("L2 distance")

pr <- ggarrange(gr1,gr2,nrow = 1,ncol = 2,labels = c("a","b"))
ggsave(pr,file="../figures/Null.pdf",width=8, height=5,scale = 0.8)


df_vuln_combined <- dplyr::inner_join(df_vuln,df_vuln_test,by="metabolite")
df_vuln_combined$diff <- df_vuln_combined$mean_abs.y - df_vuln_combined$mean_abs.x
df_vuln_combined <- df_vuln_combined[order(abs(df_vuln_combined$diff), decreasing = TRUE), ]

top_vuln <- head(df_vuln_combined, 20)


g51 <- ggplot(df_vuln, aes(x = mean_abs)) +
  geom_histogram(bins = 30, fill="#4575b4",color="white") +
  custom_theme

g52 <- ggplot(df_vuln_test, aes(x = mean_abs)) +
  geom_histogram(bins = 30, fill="#d73027",color="white") +
  custom_theme

g6 <- ggplot(top_vuln) +
  # segment (stick)
  geom_segment(data=top_vuln,aes(x = mean_abs.x, xend = mean_abs.y, y = metabolite, yend = metabolite),
    color = "grey70"
  ) +
  geom_point(data=top_vuln,aes(x=mean_abs.x,y=metabolite),color="#4575b4") +
  geom_point(data=top_vuln,aes(x=mean_abs.y,y=metabolite),color="#d73027") +
  custom_theme1 +
  xlab("Median difference (test - train)") +
  ylab("")



p2 <- ggarrange(g51,g52,nrow = 2,ncol = 1,labels = c("a","b"))
p3 <- ggarrange(p2,g6,nrow = 1,ncol = 2,labels = c("","c"),widths = c(0.5,1))

ggsave(p3,file="../figures/met_vulner.pdf",width=8.5, height=6,scale = 0.8)


g7 <- ggplot(data = meta,aes(euc,cos)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10") +
  stat_cor(method = "spearman") + custom_theme + 
  xlab("L2") + ylab("Cosine")


g8 <- ggplot(data = meta_test,aes(euc,cos)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10") +
  stat_cor(method = "spearman") + custom_theme + 
  xlab("L2") + ylab("Cosine")


p4 <- ggarrange(g7,g8,nrow = 1,ncol = 2,labels = c("a","b"))
ggsave(p4,file="../figures/L2_cosine.pdf",width=9, height=4,scale = 0.8)
