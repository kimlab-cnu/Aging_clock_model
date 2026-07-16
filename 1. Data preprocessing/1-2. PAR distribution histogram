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

# ####피규어1A. PAR 분포 히스토그램
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

###피규어1C 코드 추가
# Figure 2C. Coefficient of Variation (CV) comparison across internal standards
# 이미 계산해둔 CV 값을 막대그래프로 시각화.
# 해석 포인트: stand4, stand5 CV가 낮고 stand6 CV가 가장 큼.
cv_tbl <- read_csv("CV_Group_250605.csv") %>%
  mutate(Group = factor(Group, levels = c("PAR1_Y8","PAR2_Y6","PAR3_Y8","PAR4_Y6","PAR5_Y8","PAR6_Y6")))

ggplot(cv_tbl, aes(x = Group, y = CV, fill = Group)) +
  geom_col(width = 0.4, color = "gray40") +
  geom_text(aes(label = paste0(round(CV,1),"%")), vjust = -0.4, size = 4) +
  scale_fill_manual(values = pal) +
  labs(y = "Coefficient of Variation (%)",
       x = NULL,
       title = "CV of each internal standard") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1))

# 범례 추가
cv_tbl <- read_csv("CV_Group_250605.csv") %>%
  mutate(Group = factor(Group, 
                        levels = c("PAR1_Y8","PAR2_Y6","PAR3_Y8",
                                   "PAR4_Y6","PAR5_Y8","PAR6_Y6")))

pal <- c(
  "PAR1_Y8" = "#F4A7B9",
  "PAR2_Y6" = "#D6C26E",
  "PAR3_Y8" = "#8FD19E",
  "PAR4_Y6" = "#7DCFE0",
  "PAR5_Y8" = "#8DA9E0",
  "PAR6_Y6" = "#D79ACD"
)

ggplot(cv_tbl, aes(x = Group, y = CV, fill = Group)) +
  geom_col(width = 0.6, color = "gray40") +
  geom_text(aes(label = paste0(round(CV,1), "%")),
            vjust = -0.4, size = 4, color = "black") +
  scale_fill_manual(
    name = "Internal Standard (stand)",
    values = pal,
    labels = c("stand1 (PAR1_Y8)",
               "stand2 (PAR2_Y6)",
               "stand3 (PAR3_Y8)",
               "stand4 (PAR4_Y6)",
               "stand5 (PAR5_Y8)",
               "stand6 (PAR6_Y6)")
  ) +
  labs(
    title = "Coefficient of Variation (CV) across internal standards",
    y = "CV (%)",
    x = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",         # 범례를 오른쪽으로
    legend.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )
