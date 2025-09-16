library(limma)
library(tidyverse)
library(dplyr)

getwd()
write.csv(data_log2, "./output/Aging_PAR_log2_data.csv")

set.seed(220)


getwd()
setwd("D:/Analysis/p-value/input")
Sys.glob("*")

data <- data.frame(read.csv("Aging_PAR_log_data.csv"))

colnames(data) <- sub("^X", "Sample_", colnames(data))

row.names(data) <- data$Name
data$Name <- NULL

group <- rep(c("twenties", "thirties", "forties", "fifties"), each = 6)
metadata <- data.frame(Sample = colnames(data), Group = group)

design <- model.matrix(~ 0 + metadata$Group)
colnames(design) <- levels(factor(metadata$Group))

contrast.matrix <- makeContrasts(
  TwovsThree = thirties - twenties,
  ThreevsFour = forties - thirties,
  FourvsFive = fifties - forties,
  TwovsFour = forties - twenties,
  TwovsFive = fifties - twenties,
  ThreevsFive = fifties - thirties,
  levels = design
)                                    

data_log2 <- data * log2(10) 

fit <- lmFit(data_log2, design)

fit3 <- contrasts.fit(fit, contrast.matrix)
fit3 <- eBayes(fit3, trend=TRUE)

plotMA(fit3)
plotSA(fit3)


res_G1vsG2 <- topTable(fit3, coef = "TwovsThree", number = Inf)
res_G2vsG3 <- topTable(fit3, coef = "ThreevsFour", number = Inf)
res_G3vsG4 <- topTable(fit3, coef = "FourvsFive", number = Inf)

res_G1vsG3 <- topTable(fit3, coef = "TwovsFour", number = Inf)
res_G1vsG4 <- topTable(fit3, coef = "TwovsFive", number = Inf)
res_G2vsG4 <- topTable(fit3, coef = "ThreevsFive", number = Inf)



FC_cutoff = log2(1.5)

res_G1vsG2$Significant <- with(res_G1vsG2, ifelse(P.Value < 0.05 & abs(logFC) > FC_cutoff, 
                               ifelse(logFC > FC_cutoff, "30S", "20S"), "insig"))
res_G2vsG3$Significant <- with(res_G2vsG3, ifelse(P.Value < 0.05 & abs(logFC) > FC_cutoff, 
                               ifelse(logFC > FC_cutoff, "40S", "30S"), "insig"))
res_G3vsG4$Significant <- with(res_G3vsG4, ifelse(P.Value < 0.05 & abs(logFC) > FC_cutoff, 
                               ifelse(logFC > FC_cutoff, "50S", "40S"), "insig"))

res_G1vsG3$Significant <- with(res_G1vsG3, ifelse(P.Value < 0.05 & abs(logFC) > FC_cutoff, 
                               ifelse(logFC > FC_cutoff, "40S", "20S"), "insig"))
res_G1vsG4$Significant <- with(res_G1vsG4, ifelse(P.Value < 0.05 & abs(logFC) > FC_cutoff, 
                                                  ifelse(logFC > FC_cutoff, "50S", "20S"), "insig"))
res_G2vsG4$Significant <- with(res_G2vsG4, ifelse(P.Value < 0.05 & abs(logFC) > FC_cutoff, 
                                                  ifelse(logFC > FC_cutoff, "50S", "30S"), "insig"))

#res_G1vsG2$Significant <- with(res_G1vsG2, P.Value < 0.05 & abs(logFC) > log2(1.5))
#res_G2vsG3$Significant <- with(res_G2vsG3, P.Value < 0.05 & abs(logFC) > log2(1.5))
#res_G3vsG4$Significant <- with(res_G3vsG4, P.Value < 0.05 & abs(logFC) > log2(1.5))

ggplot(res_G1vsG2, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("blue", "red", "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: 20s vs 30s", x = "log2 Fold Change", y = "-log10(p-value)")

ggplot(res_G2vsG3, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("red", "purple", "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: 30s vs 40s", x = "log2 Fold Change", y = "-log10(p-value)")

ggplot(res_G3vsG4, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("purple", "darkgreen", "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: 40s vs 50s", x = "log2 Fold Change", y = "-log10(p-value)")

ggplot(res_G1vsG3, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("blue", "purple", "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: 20s vs 40s", x = "log2 Fold Change", y = "-log10(p-value)")

ggplot(res_G1vsG4, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("blue", "darkgreen", "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: 20s vs 50s", x = "log2 Fold Change", y = "-log10(p-value)")

ggplot(res_G2vsG4, aes(x = logFC, y = -log10(P.Value), color = Significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("red", "darkgreen", "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: 30s vs 50s", x = "log2 Fold Change", y = "-log10(p-value)")



res_G1vsG2$Ion <- rownames(res_G1vsG2)
res_G2vsG3$Ion <- rownames(res_G2vsG3)
res_G3vsG4$Ion <- rownames(res_G3vsG4)

res_G1vsG3$Ion <- rownames(res_G1vsG3)
res_G1vsG4$Ion <- rownames(res_G1vsG4)
res_G2vsG4$Ion <- rownames(res_G2vsG4)


merged_results <- res_G1vsG2 %>%
  select(Ion, logFC_20Svs30S = logFC, pval_20Svs30S = P.Value, adjp_20Svs30S = adj.P.Val) %>%
  inner_join(
    res_G2vsG3 %>%
      select(Ion, logFC_30Svs40S = logFC, pval_30Svs40S = P.Value, adjp_30Svs40S = adj.P.Val),
    by = "Ion"
  ) %>%
  inner_join(
    res_G3vsG4 %>%
      select(Ion, logFC_40Svs50S = logFC, pval_40Svs50S = P.Value, adjp_40Svs50S = adj.P.Val),
    by = "Ion"
  ) %>%
  inner_join(
    res_G1vsG3 %>%
      select(Ion, logFC_20Svs40S = logFC, pval_20Svs40S = P.Value, adjp_20Svs40S = adj.P.Val),
    by = "Ion"
  )%>%
  inner_join(
    res_G1vsG4 %>%
      select(Ion, logFC_20Svs50S = logFC, pval_20Svs50S = P.Value, adjp_20Svs50S = adj.P.Val),
    by = "Ion"
  )%>%
  inner_join(
    res_G2vsG4 %>%
      select(Ion, logFC_30Svs50S = logFC, pval_30Svs50S = P.Value, adjp_30Svs50S = adj.P.Val),
    by = "Ion"
  )



ions_uuu <- merged_results %>%
  filter(
    adjp_20Svs30S < 0.05,
    adjp_30Svs40S < 0.05,
    adjp_40Svs50S < 0.05,
    logFC_20Svs30S > FC_cutoff,   
    logFC_30Svs40S > FC_cutoff, 
    logFC_40Svs50S > FC_cutoff    
  )

ions_raw_uuu <- merged_results %>%
  filter(
    pval_20Svs30S < 0.05,
    pval_30Svs40S < 0.05,
    pval_40Svs50S < 0.05,
    logFC_20Svs30S > FC_cutoff,   
    logFC_30Svs40S > FC_cutoff, 
    logFC_40Svs50S > FC_cutoff    # 증가-증가-증가: 0개
  )

ions_ddd <- merged_results %>%
  filter(
    adjp_20Svs30S < 0.05,
    adjp_30Svs40S < 0.05,
    adjp_40Svs50S < 0.05,
    logFC_20Svs30S < -FC_cutoff,   
    logFC_30Svs40S < -FC_cutoff, 
    logFC_40Svs50S < -FC_cutoff    
  )

ions_raw_ddd <- merged_results %>%
  filter(
    pval_20Svs30S < 0.05,
    pval_30Svs40S < 0.05,
    pval_40Svs50S < 0.05,
    logFC_20Svs30S < -FC_cutoff,   
    logFC_30Svs40S < -FC_cutoff, 
    logFC_40Svs50S < -FC_cutoff    # 감소-감소-감소: 0개
  )

# pattern: 110 = (1,1,0) or (-1,-1,0)
ions_raw_uum <- merged_results %>%
  filter(
    pval_20Svs30S < 0.05,
    pval_30Svs40S < 0.05,
    pval_40Svs50S > 0.05,
    logFC_20Svs30S > FC_cutoff,   # 증가-증가-유지: 0개
    logFC_30Svs40S > FC_cutoff, 
  )

ions_raw_ddm <- merged_results %>%
  filter(
    pval_20Svs30S < 0.05,
    pval_30Svs40S < 0.05,
    pval_40Svs50S > 0.05,
    logFC_20Svs30S < -FC_cutoff,   # 감소-감소-유지: 0개
    logFC_30Svs40S < -FC_cutoff, 
  )

# pattern: 101 = (1,0,1) or (-1, 0, -1)
ions_raw_umu <- merged_results %>%
  filter(
    pval_20Svs30S < 0.05,
    pval_30Svs40S > 0.05,
    pval_40Svs50S < 0.05,
    logFC_20Svs30S > FC_cutoff,   # 증가-유지-증가 : 29개
    logFC_40Svs50S > FC_cutoff, 
  )

ions_raw_dmd <- merged_results %>%
  filter(
    pval_20Svs30S < 0.05,
    pval_30Svs40S > 0.05,
    pval_40Svs50S < 0.05,
    logFC_20Svs30S < -FC_cutoff,   # 감소-유지-감소: 0개
    logFC_40Svs50S < -FC_cutoff, 
  )

# pattern: 100 = (1,0,0) or (-1,0,0)
ions_raw_umm <- merged_results %>%
  filter(
    pval_20Svs30S < 0.05,
    pval_30Svs40S > 0.05,
    pval_40Svs50S > 0.05,
    logFC_20Svs30S > FC_cutoff,   # 증가-유지-유지 : 58개
  )

ions_raw_dmm <- merged_results %>%
  filter(
    pval_20Svs30S < 0.05,
    pval_30Svs40S > 0.05,
    pval_40Svs50S > 0.05,
    logFC_20Svs30S < -FC_cutoff,   # 감소-유지-유지 : 12개
  )

# pattern: 011 = (0,1,1) or (0,-1, -1)
ions_raw_muu <- merged_results %>%
  filter(
    pval_20Svs30S > 0.05,
    pval_30Svs40S < 0.05,
    pval_40Svs50S < 0.05,
    logFC_30Svs40S > FC_cutoff,   # 유지-증가-증가 : 2개
    logFC_40Svs50S > FC_cutoff
  )

ions_raw_mdd <- merged_results %>%
  filter(
    pval_20Svs30S > 0.05,
    pval_30Svs40S < 0.05,
    pval_40Svs50S < 0.05,
    logFC_30Svs40S < -FC_cutoff,   # 유지-감소-감소 : 0개
    logFC_40Svs50S < -FC_cutoff
  )

# pattern: 010 = (0,1,0) or (0,-1,0)
ions_raw_mum <- merged_results %>%
  filter(
    pval_20Svs30S > 0.05,
    pval_30Svs40S < 0.05,
    pval_40Svs50S > 0.05,
    logFC_30Svs40S > FC_cutoff   # 유지-증가-유지 : 71개
  )

ions_raw_mdm <- merged_results %>%
  filter(
    pval_20Svs30S > 0.05,
    pval_30Svs40S < 0.05,
    pval_40Svs50S > 0.05,
    logFC_30Svs40S < -FC_cutoff   # 유지-감소-유지 : 2개
  )

# pattern: 001 = (0,0,1) or (0,0,-1)
ions_raw_mmu <- merged_results %>%
  filter(
    pval_20Svs30S > 0.05,
    pval_30Svs40S > 0.05,
    pval_40Svs50S < 0.05,
    logFC_40Svs50S > FC_cutoff    # 유지-유지-증가 : 75개
  )

ions_raw_mmd <- merged_results %>%
  filter(
    pval_20Svs30S > 0.05,
    pval_30Svs40S > 0.05,
    pval_40Svs50S < 0.05,
    logFC_40Svs50S < -FC_cutoff   # 유지-유지-감소 : 9개
  )

# pattern: 000
ions_raw_mmm <- merged_results %>%
  filter(
    pval_20Svs30S > 0.05,
    pval_30Svs40S > 0.05,
    pval_40Svs50S > 0.05          # 유지-유지-유지 : 1521개
  )


pattern_list <- list(
  "(1,1,1)" = ions_raw_uuu,             
  "(-1,-1,-1)" = ions_raw_ddd,
  "(1,1,0)" = ions_raw_uum,
  "(-1,-1,0)" = ions_raw_ddm,
  "(1,0,1)" = ions_raw_umu,
  "(-1,0,-1)" = ions_raw_dmd,#
  "(1,0,0)" = ions_raw_umm,
  "(-1,0,0)" = ions_raw_dmm,#
  "(0,1,1)" = ions_raw_muu,
  "(0,-1,-1)" = ions_raw_mdd,
  "(0,1,0)" = ions_raw_mum,
  "(0,-1,0)" = ions_raw_mdm,
  "(0,0,1)" = ions_raw_mmu,
  "(0,0,-1)" = ions_raw_mmd,
  "(0,0,0)" = ions_raw_mmm
)


pattern_ions <- purrr::imap_dfr(pattern_list, ~ mutate(.x, pattern = .y))  

getwd()
setwd("./output/")

write.csv(pattern_ions, "ion_pattern_p-val_and_logFC_fit3_limma.csv", row.names=FALSE)


ggplot(pattern_ions, aes(x = pattern)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  labs(title = "Number of ions per pattern", y = "Ion count", x = "Pattern") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


pattern_ions_sig <- pattern_ions %>%
  filter(pattern != "(0,0,0)")


data_long <- data_log2 %>%
  rownames_to_column("Ion") %>%
  pivot_longer(-Ion, names_to = "Sample", values_to = "PAR")


data_long <- data_long %>%
  left_join(metadata, by = "Sample")

data_pattern <- data_long %>%
  inner_join(pattern_ions_sig %>% select(Ion, pattern), by = "Ion")


mean_trends <- data_pattern %>%
  group_by(pattern, Group, Ion) %>%
  summarise(PAR_mean = mean(PAR, na.rm = TRUE), .groups = "drop")

mean_trends$Group <- factor(mean_trends$Group, levels = c("twenties", "thirties", "forties", "fifties"))

ggplot(mean_trends, aes(x = Group, y = PAR_mean, group = Ion)) +
  geom_line(alpha = 0.3) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", linewidth = 1.2) +
  facet_wrap(~ pattern, scales = "free_y") +
  theme_minimal() +
  labs(title = "Groupwise Expression Trajectories per Pattern", y = "Mean PAR")




library(pheatmap)

sig_2030s_ions <- res_G1vsG2 %>% filter(Significant != "insig") %>% pull(Ion)
sig_2040s_ions <- res_G1vsG3 %>% filter(Significant != "insig") %>% pull(Ion)
sig_2050s_ions <- res_G1vsG4 %>% filter(Significant != "insig") %>% pull(Ion)
sig_3050s_ions <- res_G2vsG4 %>% filter(Significant != "insig") %>% pull(Ion)
sample_2030s <- metadata %>% filter(Group %in% c("twenties", "thirties")) %>% pull(Sample)

heatmap_2030s <- data_log2[sig_2030s_ions, sample_2030s]
pheatmap(heatmap_2030s, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames= FALSE, scale="row")

sig_3040s_ions <- res_G2vsG3 %>% filter(Significant != "insig") %>% pull(Ion)
sample_3040s <- metadata %>% filter(Group %in% c("thirties", "forties")) %>% pull(Sample)

heatmap_3040s <- data_log2[sig_3040s_ions, sample_3040s]
pheatmap(heatmap_3040s, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames= FALSE, scale="row")

sig_4050s_ions <- res_G3vsG4 %>% filter(Significant != "insig") %>% pull(Ion)
sample_4050s <- metadata %>% filter(Group %in% c("forties", "fifties")) %>% pull(Sample)

heatmap_4050s <- data_log2[sig_4050s_ions, sample_4050s]
pheatmap(heatmap_4050s, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames= FALSE, scale="row")

library(gridExtra)
library(grid)   # 제목 넣기 

p1 <- pheatmap(heatmap_2030s, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames= FALSE, legend = FALSE, scale="row")
p2 <- pheatmap(heatmap_3040s, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames= FALSE, legend = FALSE, scale="row")
p3 <- pheatmap(heatmap_4050s, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames= FALSE, scale="row")

grid.arrange(grobs = list(p1[[4]], p2[[4]], p3[[4]]), ncol = 3, 
             top = textGrob("Differentially Expressed Ions", gp = gpar(fontsize = 16, fontface = "bold")))

all_ions <- unique(c(sig_2030s_ions, sig_2040s_ions, sig_2050s_ions))
venn_data <- data.frame(
  Comp_20Svs30S = all_ions %in% sig_2030s_ions,
  Comp_20Svs40S = all_ions %in% sig_2040s_ions,
  Comp_20Svs50S = all_ions %in% sig_2050s_ions
)

rownames(venn_data) <- all_ions

vc <- vennCounts(venn_data)
vennDiagram(vc, main = "DE Ion Overlap")

common_ions <- Reduce(intersect, list(sig_2030s_ions, sig_3040s_ions, sig_4050s_ions))

df_common <- data_long %>% filter(Ion %in% common_ions)
df_common$Group <- factor(df_common$Group, levels = c("twenties", "thirties", "forties", "fifties"))

library(ggpubr)

ggplot(df_common, aes(x = Group, y = PAR, fill = Group)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  stat_compare_means(method = "wilcox.test",  
                     comparisons = list(c("twenties", "thirties"),
                                        c("thirties", "forties"), c("forties", "fifties")), 
                    label = "p.signif") +
  facet_wrap(~ Ion, scales = "free_y") +
  theme_minimal() +
  labs(title = "Boxplot of Common DE Ions", y = "Expression (log2)", x = "")

