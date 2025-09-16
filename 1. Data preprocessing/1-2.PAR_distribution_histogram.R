getwd()

setwd("C:/Users/kjkh0/Desktop/Analysis/PAR_calculation/output/")

Sys.glob("*")

df = data.frame(read.csv("Transition_PAR_Results_250605.csv"))

head(df, 5)

PAR_data = data.frame(
  PAR = c(df$PAR1, df$PAR2, df$PAR3, df$PAR4, df$PAR5, df$PAR6),
  Group = factor(rep(c("PAR1_Y8", "PAR2_Y6", "PAR3_Y8", "PAR4_Y6", "PAR5_Y8", "PAR6_Y6"), each=45096))
)

PAR_data$PAR <- log10(PAR_data$PAR)                
PAR_data <- PAR_data[is.finite(PAR_data$PAR), ]    

df2 = data.frame(read.csv("CV_Group_250605.csv"))

df2_CV <- df2 %>%
  mutate(CV_label = paste0(Group, " (", round(CV, 1), "%)")) %>%
  arrange(Group)

PAR_data <- PAR_data %>%
  left_join(df2_CV %>% select(Group, CV_label), by = "Group")

library(ggplot2)

p <- ggplot(PAR_data, aes(x = PAR, fill = CV_label)) +
  annotate("rect", xmin=-1, xmax=1, ymin=0, ymax=Inf, 
           alpha=0.15, fill="yellow") +
  geom_histogram(position = "identity", binwidth=0.1, alpha=0.5, color="black") +
  facet_wrap(~ Group, scale="free_y")+
  coord_cartesian(xlim = c(-2, 2)) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(
    title = "Peak Area Ratio Distribution",
    x = "PAR", y = "Frequency",
    fill = "Group (CV %)"
  )

p

library(dplyr)

PAR_data %>% filter(PAR >= 0.1 & PAR <= 10) %>%
  group_by(Group) %>% summarise(count = n()) %>%
  arrange(desc(count))

PAR_data %>% filter(PAR >= 0 & PAR < 0.1) %>%
  group_by(Group) %>% summarise(count = n()) %>%
  arrange(desc(count))

PAR_data %>% filter(PAR > 10) %>%
  group_by(Group) %>% summarise(count = n()) %>%
  arrange(desc(count))
