library(dplyr)
library(tidyr)

getwd()

setwd("C:/Users/USER/Desktop/Analysis/p-value/input/")

Sys.glob("*")

data <- data.frame(read.csv("Aging_PAR_log_norm_data.csv"))

head(data, 3)

rownames(data) <- data$Name

colnames(data) <- sub("^X", "Sample_", colnames(data))

data_long <- data %>%
  pivot_longer(-Name, names_to = "sample", values_to = "value") %>%
  mutate(group = gsub("\\.\\d+$", "", sample))  

data_avg <- data_long %>%
  group_by(Name, group) %>%      
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = value) %>%
  tibble::column_to_rownames("Name")

#write.csv(data_avg, file = "PAR_log_norm_avg.csv")



library(tidyverse)

## 예상 패턴 벡터 (x, y, z)
# 감소: -1, 유지: 0, 증가: 1
# 3*3*3=27가지
# 1로 시작: 1.(1, −1, −1), 2.(1, −1, 0), 3.(1, −1, 1), 4.(1, 0, −1), 5.(1, 0, 0), 6.(1, 0, 1), 7.(1, 1, −1), 8.(1, 1, 0)*, 9.(1, 1, 1)*
# 0으로 시작: 1.(0, −1, −1), 2.(0, −1, 0), 3.(0, −1, 1), 4.(0, 0, −1), 5.(0, 0, 0), 6.(0, 0, 1), 7.(0, 1, −1), 8.(0, 1, 0), 9.(0, 1, 1)
# -1로 시작: 1.(−1, −1, −1)*, 2.(−1, −1, 0)*, 3.(−1, −1, 1), 4.(−1, 0, −1), 5.(−1, 0, 0), 6.(−1, 0, 1), 7.(−1, 1, −1), 8.(−1, 1, 0), 9.(−1, 1, 1)

df_wide <- read.csv("PAR_log_norm_avg.csv") 

df_wide <- df_wide %>%
  rename(Name = X)

target_vec <- c(1, 1, 1)


df_cosine <- df_wide %>%
  rowwise() %>%
  mutate(
    vec = list(c(`Sample_30s` - `Sample_20s`, `Sample_40s` - `Sample_30s`, `Sample_50s` - `Sample_40s`)),
    cosine = sum(unlist(vec) * target_vec) /
      (sqrt(sum(unlist(vec)^2)) * sqrt(sum(target_vec^2)))
  ) %>%
  ungroup() %>%
  select(Name, cosine) %>%
  arrange(desc(cosine))

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
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


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



top_features <- df_cosine %>%
  filter(cosine >= 0.9) %>%
  pull(Name)

df_top_long <- df_wide %>%
  filter(Name %in% top_features) %>%
  pivot_longer(cols = c(`Sample_20s`, `Sample_30s`, `Sample_40s`, `Sample_50s`),
               names_to = "group", values_to = "par") %>%
  mutate(group = factor(group, levels = c("Sample_20s", "Sample_30s", "Sample_40s", "Sample_50s")))


ggplot(df_top_long, aes(x = group, y = par, group = Name)) +
  geom_line(alpha = 0.6, color = "steelblue") +
  labs(title = "Expression Trajectories of Features with Cosine ≥ 0.9",
       x = "Group", y = "Mean Peak Area Ratio(PAR)") +
  theme_minimal()


library(pheatmap)


df_mat <- df_wide %>%
  filter(Name %in% top_features) %>%
  column_to_rownames("Name") %>%
  as.matrix()

pheatmap(df_mat,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         main = "Cluster 1 – Cosine ≥ 0.9 Features")


#write.csv(df, "output.csv", row.names = FALSE, fileEncoding = "UTF-8")

df_top_only <- df_wide %>% 
  filter(Name %in% top_features)

write.csv(df_top_only, "cosine_0.9_targets_(-1,1,1).csv", row.names = FALSE)
