# data preprocessing
library(dplyr)
library(tidyr)
library(readr)

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

# 개별샘플을 Sample_20s, 30s, 40s, 50s로 묶어 평균 PAR 구함 (코사인 유사도 볼때는 평균값이 필요)
data_avg <- data_long %>%
  group_by(Name, group) %>%      
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%   # 그룹별 평균을 구함 
  pivot_wider(names_from = group, values_from = value) %>%             # Wide 형태로 변환
  tibble::column_to_rownames("Name")  # 다시 rownames로 복원

## fold_change 용 raw data 만들기
setwd("Z:/Project/Aging_clock_model/Data/PAR_p-value/input")
Sys.glob("*.csv")

data <- data.frame(read.csv("PAR_raw_data.csv"))  

data_raw <- data %>%
  group_by(Name, Replicate) %>%
  summarise(PAR = mean(PAR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Replicate, values_from = PAR)

#data_raw$Name <- NULL

write_csv(data_raw, "./Aging_PAR_raw_data.csv")

# TITIN 중 세 이온이 결과가 2개씩 존재함을 확인 (par값은 같으므로, 평균내서 하나로만 사용하기. 아마 스카이라인 프로그램 문제인 듯)
a <- data %>%
  group_by(Name, Replicate) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)
