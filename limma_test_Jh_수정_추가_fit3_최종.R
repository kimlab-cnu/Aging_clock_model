리드미는 한국어로
코드 주석은 영어로?
주석은 지우거나 영어로

library(limma)
library(tidyverse)
library(dplyr)

# ------------- R 이미지 저장 (변수, 패키지 등 환경을 빠르게 불러올 수 있음) -------------
#save.image("../rdata/Aging_Limma_test.RData")
#load("C:/Users/USER/Desktop/Analysis/PAR_p-value/rdata/Aging_Limma_test.RData")
# ----------------------------------------------------------------------------------------

getwd()
write.csv(data_log2, "./output/Aging_PAR_log2_data.csv")

# 통계 분석을 위한 Seed (재현성 확보) 
set.seed(220)

#setwd("../input/")  # Analysis > PAR_p-value > input 경로로 들어가도록 설정해주기 

getwd()
#setwd("../")
setwd("D:/Analysis/p-value/input")
Sys.glob("*")

data <- data.frame(read.csv("Aging_PAR_log_data.csv"))

# R에서 열, 변수 이름은 숫자로 시작하면 안돼서 x가 자동으로 들어감 -> 이를 Sample_20s.1로 변경
colnames(data) <- sub("^X", "Sample_", colnames(data))

# Ion 이름은 문자열(str)이니 분석을 위해서 row.names로 빼주고 데이터는 날림림
row.names(data) <- data$Name
data$Name <- NULL

# Sample의 그룹 설정, 이를 통해 metadata를 만들어 비교군을 연결해줄 수 있음 
group <- rep(c("twenties", "thirties", "forties", "fifties"), each = 6)
metadata <- data.frame(Sample = colnames(data), Group = group)

# metadata를 활용해서 실험 디자인을 분석에 녹여낼 수 있음 
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
)                                    # 이렇게 바꿔야 20대-30대, 30대-40대, 40대-50대 비교가 됨.

# limma는 log2 transformed data에 잘 작동 -> log10 transform 결과를 log2로 변환
data_log2 <- data * log2(10) 

# --------------------------------------------------------------------

fit <- lmFit(data_log2, design)                  # 선형 모델 적합
#fit2 <- contrasts.fit(fit, contrast.matrix)  # 대비(contrast) 적용
#fit2 <- eBayes(fit2)                        # Empirical Bayes 보정

# -------------------------- Test ----------------------------------
fit3 <- contrasts.fit(fit, contrast.matrix)  # 대비(contrast) 적용
fit3 <- eBayes(fit3, trend=TRUE)            # Empirical Bayes 보정

# trend=TRUE: limma-trend; intensity 의존적 variance가 있는 MS 데이터에 적합.
# intensity에 따른 분산 변화(trend)를 반영해주기에 더 좋은 모델링이 가능하게 됨. 

plotMA(fit2)  # 평균 표현량(A), log2 Fold Change(M)으로 구성; 차등 발현 확인 및 시각화
plotSA(fit2)  # 평균 표현랑(A), 잔차 표준편차(SD)로 구성; 모델 품질과 정규화 상태 진단

plotMA(fit3)
plotSA(fit3)

# ------------------------------ 결과 해석 -------------------------------------
# fit2는 trend=TRUE 미적용 결과 -> 파란선이 낮게 깔리며 파란선 위의 점은 high noise로 간주됨
# fit3는 trend=TRUE 적용 결과 -> 파란선이 그래도 expression pattern을 따라 가는 형태 
# 그러나, 0쪽 구간에 점 들이 많이 몰려 있어 low-intensity (high-noise)가 많은 것으로 보임
# ------------------------------------------------------------------------------

# Group1 vs Group2 결과 예시
res_G1vsG2 <- topTable(fit3, coef = "TwovsThree", number = Inf)
res_G2vsG3 <- topTable(fit3, coef = "ThreevsFour", number = Inf)
res_G3vsG4 <- topTable(fit3, coef = "FourvsFive", number = Inf)

res_G1vsG3 <- topTable(fit3, coef = "TwovsFour", number = Inf)
res_G1vsG4 <- topTable(fit3, coef = "TwovsFive", number = Inf)
res_G2vsG4 <- topTable(fit3, coef = "ThreevsFive", number = Inf)

# -----------------------------------------------------------------------------------------
# adjusted p-value 기법들 확인하기 
#res_G1vsG2$adj_BH <- p.adjust(res_G1vsG2$P.Value, method = "BH")
#res_G1vsG2$adj_BF <- p.adjust(res_G1vsG2$P.Value, method = "bonferroni")   # 보수적인 조정 방법
#res_G1vsG2$adj_BY <- p.adjust(res_G1vsG2$P.Value, method = "BY")           # Benjamini-Yekutieli (더 보수적)
#res_G1vsG2$adj_holm <- p.adjust(res_G1vsG2$P.Value, method = "holm")         # holm's method (BF보다 덜 보수적)
#res_G1vsG2$fdr <- p.adjust(res_G1vsG2$P.Value, method="fdr")
# default로 적용되는 BH (=FDR)이 가장 rough한 적용 방법임에도 adjusted p-value가 높음 
# -----------------------------------------------------------------------------------------

FC_cutoff = log2(1.5)
# volcano plot 그리기
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
  scale_color_manual(values = c("red", "purple", "grey")) +   # 2개 일때는 "grey", "red" (insig, sig)
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


# 해당 데이터 프레임의 Ion 이름들을 rownames로 변환 
res_G1vsG2$Ion <- rownames(res_G1vsG2)
res_G2vsG3$Ion <- rownames(res_G2vsG3)
res_G3vsG4$Ion <- rownames(res_G3vsG4)

res_G1vsG3$Ion <- rownames(res_G1vsG3)
res_G1vsG4$Ion <- rownames(res_G1vsG4)
res_G2vsG4$Ion <- rownames(res_G2vsG4)

# 필요한 Ion, logFC, P.Value, adj.P.Val만 취해서 하나의 데이터 프레임으로 합쳐줌 
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

# ------------------------ 패턴별 이온 확인 ---------------------------
# 패턴 정의 (u: up (증가), m: maintenance (유지), d: down (감소))
# ex. mum : 유지-증가-유지
# 필터링 조건은 p-value와 logFC를 동시에 사용 (유지일 때는 p-value만 사용, logFC 조건이 딱히 없음)

# pattern: 111 = (1,1,1) or (-1,-1,-1)
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
    logFC_20Svs30S < -FC_cutoff,   # 감소-감소-유지: 0개 ###
    logFC_30Svs40S < -FC_cutoff, 
  )

# pattern: 101 = (1,0,1) or (-1, 0, -1)
ions_raw_umu <- merged_results %>%
  filter(
    pval_20Svs30S < 0.05,
    pval_30Svs40S > 0.05,
    pval_40Svs50S < 0.05,
    logFC_20Svs30S > FC_cutoff,   # 증가-유지-증가 : 31개(fit2) / 29개(fit3)
    logFC_40Svs50S > FC_cutoff, 
  )

ions_raw_dmd <- merged_results %>%
  filter(
    pval_20Svs30S < 0.05,
    pval_30Svs40S > 0.05,
    pval_40Svs50S < 0.05,
    logFC_20Svs30S < -FC_cutoff,   # 감소-유지-감소: 0개(fit2, fit3)
    logFC_40Svs50S < -FC_cutoff, 
  )

# pattern: 100 = (1,0,0) or (-1,0,0)
ions_raw_umm <- merged_results %>%
  filter(
    pval_20Svs30S < 0.05,
    pval_30Svs40S > 0.05,
    pval_40Svs50S > 0.05,
    logFC_20Svs30S > FC_cutoff,   # 증가-유지-유지 : 57개 (fit2)/ 58개 (fit3)
  )

ions_raw_dmm <- merged_results %>%
  filter(
    pval_20Svs30S < 0.05,
    pval_30Svs40S > 0.05,
    pval_40Svs50S > 0.05,
    logFC_20Svs30S < -FC_cutoff,   # 감소-유지-유지 : 11개 (fit2)/ 12개 (fit3)
  )

# pattern: 011 = (0,1,1) or (0,-1, -1)
ions_raw_muu <- merged_results %>%
  filter(
    pval_20Svs30S > 0.05,
    pval_30Svs40S < 0.05,
    pval_40Svs50S < 0.05,
    logFC_30Svs40S > FC_cutoff,   # 유지-증가-증가 : 2개 (fit2, fit3)
    logFC_40Svs50S > FC_cutoff
  )

ions_raw_mdd <- merged_results %>%
  filter(
    pval_20Svs30S > 0.05,
    pval_30Svs40S < 0.05,
    pval_40Svs50S < 0.05,
    logFC_30Svs40S < -FC_cutoff,   # 유지-감소-감소 : 0개 (fit2, fit3)
    logFC_40Svs50S < -FC_cutoff
  )

# pattern: 010 = (0,1,0) or (0,-1,0)
ions_raw_mum <- merged_results %>%
  filter(
    pval_20Svs30S > 0.05,
    pval_30Svs40S < 0.05,
    pval_40Svs50S > 0.05,
    logFC_30Svs40S > FC_cutoff   # 유지-증가-유지 : 74개 (fit2) / 71개 (fit3)
  )

ions_raw_mdm <- merged_results %>%
  filter(
    pval_20Svs30S > 0.05,
    pval_30Svs40S < 0.05,
    pval_40Svs50S > 0.05,
    logFC_30Svs40S < -FC_cutoff   # 유지-감소-유지 : 2개 (fit2, fit3)
  )

# pattern: 001 = (0,0,1) or (0,0,-1)
ions_raw_mmu <- merged_results %>%
  filter(
    pval_20Svs30S > 0.05,
    pval_30Svs40S > 0.05,
    pval_40Svs50S < 0.05,
    logFC_40Svs50S > FC_cutoff    # 유지-유지-증가 : 79개 (fit2) / 75개 (fit3)
  )

ions_raw_mmd <- merged_results %>%
  filter(
    pval_20Svs30S > 0.05,
    pval_30Svs40S > 0.05,
    pval_40Svs50S < 0.05,
    logFC_40Svs50S < -FC_cutoff   # 유지-유지-감소 : 9개 (fit2, fit3)
  )

# pattern: 000
ions_raw_mmm <- merged_results %>%
  filter(
    pval_20Svs30S > 0.05,
    pval_30Svs40S > 0.05,
    pval_40Svs50S > 0.05          # 유지-유지-유지 : 1545개 (fit2) / 1521개 (fit3)
  )

##### 전체 pattern 리스트를 넣어줌 #####
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

# pattern 리스트를 정의하고, 각 데이터별로 pattern 열을 추가해서 패턴 정보를 넣어줌.
pattern_ions <- purrr::imap_dfr(pattern_list, ~ mutate(.x, pattern = .y))  
# 65개는 어디서 놓친걸까? -> 유의미하지 않은 패턴 (W, M 모양)으로 보임 


#intersect(cosine_data$Name, ion_foldchange$Ion)


getwd()
setwd("./output/")

write.csv(pattern_ions, "ion_pattern_p-val_and_logFC_fit3_limma.csv", row.names=FALSE)

## downstream 시각화
ggplot(pattern_ions, aes(x = pattern)) +
  geom_bar(fill = "steelblue") +
  theme_minimal() +
  labs(title = "Number of ions per pattern", y = "Ion count", x = "Pattern") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 활용 예시
# (0,0,0) = 유지/유지/유지가 아닌 유의미한 패턴을 가지는 이온만 선별
pattern_ions_sig <- pattern_ions %>%
  filter(pattern != "(0,0,0)")

# 다시 long form 형태의 dataframe을 구성
data_long <- data_log2 %>%
  rownames_to_column("Ion") %>%
  pivot_longer(-Ion, names_to = "Sample", values_to = "PAR")

# metadata에서 Sample 정보를 통해 Group 정보를 연결해줌 (같은 key -> join 가능)
data_long <- data_long %>%
  left_join(metadata, by = "Sample")

data_pattern <- data_long %>%
  inner_join(pattern_ions_sig %>% select(Ion, pattern), by = "Ion")

# PAR 평균을 구해서 패턴을 그림 (K-means 때와 동일)
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

# 위에서 정의한 총 9개 패턴에 대한 시각화 결과. ion 개수도 나오고, 평균 선도 그려줌. 

# ----------------------- downstream 시각화 (figure 그리기) ---------------------------
# Volcano plot → DE protein 확인
# Barplot → Up/Downregulated 수 요약
# Heatmap → DE ion intensity 패턴
# Boxplot → 대표 ion/group 분포 비교
# Pathway enrichment → dotplot/emapplot
# PCA → 전체 샘플/DE ion 기준 군집 확인
# Cosine trend line plot → 패턴 분리된 ion 표현
# ridge/violin plot → 각 그룹의 분포 감각적으로 보여줌

#install.packages("pheatmap")     # 패키지 없는 경우 - # 제거 하고 설치
library(pheatmap)

# ----- 방식 -----
# 1. signifcant한 iont만을 추출함 (filter로 insig가 아닌 것들을 뽑음)
# 2. metadata에서 해당 비교군에 맞는 sample 정보만 가져옴
# 3. heatmap을 그리기 위한 데이터를 data_log2 (limma 입력 데이터)에서 추출함
# 4. pheatmap 라이브러리로 heatmap을 그림
# ----------------

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

#install.packages("gridExtra")
library(gridExtra)
library(grid)   # 제목 넣기 

## heatmap 여러 개를 한번에 배치하는 방법 (다른 plot에도 동일하게 적용 가능)
p1 <- pheatmap(heatmap_2030s, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames= FALSE, legend = FALSE, scale="row")
p2 <- pheatmap(heatmap_3040s, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames= FALSE, legend = FALSE, scale="row")
p3 <- pheatmap(heatmap_4050s, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames= FALSE, scale="row")

grid.arrange(grobs = list(p1[[4]], p2[[4]], p3[[4]]), ncol = 3, 
             top = textGrob("Differentially Expressed Ions", gp = gpar(fontsize = 16, fontface = "bold")))

## venn-diagram (겹치는 ion이 있는지)
all_ions <- unique(c(sig_2030s_ions, sig_2040s_ions, sig_2050s_ions))
venn_data <- data.frame(
  Comp_20Svs30S = all_ions %in% sig_2030s_ions,
  Comp_20Svs40S = all_ions %in% sig_2040s_ions,
  Comp_20Svs50S = all_ions %in% sig_2050s_ions
)

rownames(venn_data) <- all_ions

vc <- vennCounts(venn_data)
vennDiagram(vc, main = "DE Ion Overlap")

# ------------- venn diagram을 어떻게 해석할 것인가? ------------------------------- #
common_ions <- Reduce(intersect, list(sig_2030s_ions, sig_3040s_ions, sig_4050s_ions))

df_common <- data_long %>% filter(Ion %in% common_ions)
df_common$Group <- factor(df_common$Group, levels = c("twenties", "thirties", "forties", "fifties"))

## boxplot (개별 이온 plotting)
library(ggpubr)  # 통계분석 추가 (p-value)

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


