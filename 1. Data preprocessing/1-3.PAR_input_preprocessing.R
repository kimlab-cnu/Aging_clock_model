getwd()

setwd("../../Desktop/Analysis/PAR_Kmeans_clustering/input/")

Sys.glob("*")

data <- data.frame(read.csv('PAR_raw_data.csv'))

library(dplyr)

data <- data %>% arrange(Name)

sample_list <- factor(data$Replicate)
group_list <- data$Group

library(tidyr)

data2 <- data %>% pivot_wider(id_cols = c(Name), names_from = Replicate, values_from = PAR,
                              values_fn = list(PAR = mean)) 
  
write.csv(data2, "./PAR_Kmeans_input.csv")

group_list <- c("group", rep("20S", 6), rep("30S", 6), rep("40S", 6), rep("50S", 6))

colnames(data2)

data3 <- rbind(group_list, data2)

write.csv(data3, "./PAR_Kmeans_group.csv")