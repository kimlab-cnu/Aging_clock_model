# data preprocessing
library(dplyr)
library(tidyr)

getwd()

setwd("C:/Users/USER/Desktop/Analysis/p-value/input/")

Sys.glob("*")

#data <- read_csv("Aging_PAR_log_norm_data.csv", show_col_types=FALSE) 
# csv로 읽어서 열에 숫자 이름(20s-1)이 유지되도록 하는 방법
data <- data.frame(read.csv("Aging_PAR_log_norm_data.csv"))
# log transformation + z-score normalization 데이터 사용

head(data, 3)

## data.frame
rownames(data) <- data$Name

colnames(data) <- sub("^X", "Sample_", colnames(data))

data_long <- data %>%
  pivot_longer(-Name, names_to = "sample", values_to = "value") %>%
    # 모든 열을 sample로, 각 sample별 값은 value로 저장. (새로운 열을 만듬)
  mutate(group = gsub("\\.\\d+$", "", sample))  
    # group 열에 Sample_20s.1 -> Sample_20s로 변환해서 추가 (group 정보)

# 개별샘플을 Sample_20s, 30s, 40s, 50s로 묶어 평균 PAR 구함
data_avg <- data_long %>%
  group_by(Name, group) %>%      
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%   # 그룹별 평균을 구함 
  pivot_wider(names_from = group, values_from = value) %>%             # Wide 형태로 변환
  tibble::column_to_rownames("Name")  # 다시 rownames로 복원

#write.csv(data_avg, file = "PAR_log_norm_avg.csv")

##########################################
# 그룹별 평균 PAR 데이터로 코사인 유사도 #
##########################################
library(tidyverse)

## 예상 패턴 벡터 (x, y, z)
# 감소: -1, 유지: 0, 증가: 1
# 3*3*3=27가지
# 1로 시작: 1.(1, −1, −1), 2.(1, −1, 0), 3.(1, −1, 1), 4.(1, 0, −1), 5.(1, 0, 0), 6.(1, 0, 1), 7.(1, 1, −1), 8.(1, 1, 0)*, 9.(1, 1, 1)*
# 0으로 시작: 1.(0, −1, −1), 2.(0, −1, 0), 3.(0, −1, 1), 4.(0, 0, −1), 5.(0, 0, 0), 6.(0, 0, 1), 7.(0, 1, −1), 8.(0, 1, 0), 9.(0, 1, 1)
# -1로 시작: 1.(−1, −1, −1)*, 2.(−1, −1, 0)*, 3.(−1, −1, 1), 4.(−1, 0, −1), 5.(−1, 0, 0), 6.(−1, 0, 1), 7.(−1, 1, −1), 8.(−1, 1, 0), 9.(−1, 1, 1)

# 전체 평균값 데이터 불러오기 
df_wide <- read.csv("PAR_log_norm_avg.csv") 

df_wide <- df_wide %>%                 # 기존 데이터프레임
  rename(Name = X)           # X → Name

target_vec <- c(1, 1, 1)

#data_avg

# 각 feature별 실제 변화 벡터 계산 & cosine 유사도 계산
df_cosine <- df_wide %>%
  rowwise() %>%
  mutate(
    vec = list(c(`Sample_30s` - `Sample_20s`, `Sample_40s` - `Sample_30s`, `Sample_50s` - `Sample_40s`)),
    cosine = sum(unlist(vec) * target_vec) /
      (sqrt(sum(unlist(vec)^2)) * sqrt(sum(target_vec^2)))
  ) %>%
  ungroup() %>%
  select(Name, cosine) %>%
  arrange(desc(cosine))  # 유사도 높은 순 정렬

library(ggplot2)

ggplot(df_cosine, aes(x = reorder(Name, -cosine), y = cosine)) +
  geom_point(aes(color = cosine >= 0.9), size = 3) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "red", linewidth = 1) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "blue")) +
  labs(title = "Cosine Similarity to Expected Pattern",
       x = "Feature (Peptide Name)",
       y = "Cosine Similarity",
       color = "≥ 0.9") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # 많으면 생략 가능
        axis.ticks.x = element_blank())

# cosine similarity ≥ 0.9인 feature만 추출
df_cosine_filtered <- df_cosine %>%
  filter(cosine >= 0.9) %>%
  arrange(desc(cosine))

print(df_cosine_filtered, n = Inf)

ggplot(df_cosine, aes(x = cosine)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "white") +
  geom_vline(xintercept = 0.9, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Cosine Similarity to Target Pattern",
       x = "Cosine Similarity", y = "Count")

mean(df_cosine$cosine >= 0.9)
mean(df_cosine$cosine >= 0.85)
mean(df_cosine$cosine >= 0.8)
mean(df_cosine$cosine >= 0.7)

#기준점	해당 비율	해석 (K0~)
#≥ 0.9	3.7% |23.8% | 
#≥ 0.85	2.4% |8.5% | 
#≥ 0.8	3.1% |10.9%| 
#≥ 0.7	6.2% |13.4%|
#

### cosine 0.8 이상 feature를 추출
## cosine similarity : 0이면 orthogonal (무관), 1이면 유사도 100%
# 보통은 0.7이상 유사, 0.8이상 상당히 유사, 0.9이상 엄격하게 적용
# histogram과 mean 함수를 통해 cosine similarity 분포를 확인하고 진행
top_features <- df_cosine %>%
  filter(cosine >= 0.9) %>%
  pull(Name)

df_top_long <- df_wide %>%
  filter(Name %in% top_features) %>%
  pivot_longer(cols = c(`Sample_20s`, `Sample_30s`, `Sample_40s`, `Sample_50s`),
               names_to = "group", values_to = "par") %>%
  mutate(group = factor(group, levels = c("Sample_20s", "Sample_30s", "Sample_40s", "Sample_50s")))

# lineplot으로 시각화# K0_lineplot_Expression Trajectories of Features with Cosine ≥ 0.9
ggplot(df_top_long, aes(x = group, y = par, group = Name)) +
  geom_line(alpha = 0.6, color = "steelblue") +
  labs(title = "Expression Trajectories of Features with Cosine ≥ 0.9",
       x = "Group", y = "Mean Peak Area Ratio(PAR)") +
  theme_minimal()

##### heatmap 형태로 패턴 양상 확인### K0_heatmap_Cluster 1 – Cosine ≥ 0.9 Features
library(pheatmap)

# wide → matrix 형태로 변환
df_mat <- df_wide %>%
  filter(Name %in% top_features) %>%
  column_to_rownames("Name") %>%
  as.matrix()

pheatmap(df_mat,
         scale = "row",  # z-score 정규화로 패턴 강조
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         main = "Cluster 1 – Cosine ≥ 0.9 Features")


#write.csv(df, "output.csv", row.names = FALSE, fileEncoding = "UTF-8")



# wide 데이터에서 해당 타깃 행만 추출
df_top_only <- df_wide %>% 
  filter(Name %in% top_features)

# 필요하면 파일로 저장
write.csv(df_top_only, "cosine_0.9_targets_(-1,1,1).csv", row.names = FALSE)


# 전체 데이터에서 p-value 구하고 그 결과를 코사인 유사도와 비교 (겹치는게 있는지)

