library(dplyr)
library(tidyr)
library(readr)

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


setwd("Z:/Project/Aging_clock_model/Data/PAR_p-value/input")
Sys.glob("*.csv")

data <- data.frame(read.csv("PAR_raw_data.csv"))  

data_raw <- data %>%
  group_by(Name, Replicate) %>%
  summarise(PAR = mean(PAR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Replicate, values_from = PAR)

write_csv(data_raw, "./Aging_PAR_raw_data.csv")



a <- data %>%
  group_by(Name, Replicate) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)
