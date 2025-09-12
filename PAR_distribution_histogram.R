getwd()

setwd("C:/Users/kjkh0/Desktop/Analysis/PAR_calculation/output/")

Sys.glob("*")  # 해당 경로에 어떤 파일들이 있는지 확인하는 방법

#########
df = data.frame(read.csv("Transition_PAR_Results_250605.csv"))
# Data_preprocessing_for_PAR.R 파일을 통해 수행한 결과 (PAR1, 2, 3)

head(df, 5)  # 열 정보와 데이터가 잘 기입되었는지, NA가 있는지 없는지 확인 가능

PAR_data = data.frame(
  PAR = c(df$PAR1, df$PAR2, df$PAR3, df$PAR4, df$PAR5, df$PAR6),
  Group = factor(rep(c("PAR1_Y8", "PAR2_Y6", "PAR3_Y8", "PAR4_Y6", "PAR5_Y8", "PAR6_Y6"), each=45096))
)

PAR_data$PAR <- log10(PAR_data$PAR)                # 업데이트 하기 전에는 변수 이름 바꿔서 확인해보기기
PAR_data <- PAR_data[is.finite(PAR_data$PAR), ]    # -Inf가 너무 많아서 잘라주는 역할

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

# PAR1_Y8 CV = ~~ %
# PAR1_Y8: ~~ %

library(dplyr)

PAR_data %>% filter(PAR >= 0.1 & PAR <= 10) %>%
  group_by(Group) %>% summarise(count = n()) %>%
  arrange(desc(count))

PAR_data %>% filter(PAR >= 0 & PAR < 0.1) %>%
  group_by(Group) %>% summarise(count = n()) %>%
  arrange(desc(count))

##### A tibble : 4 x 3########################
## Group    0.1 <= x <= 10  |  0 <= x < 0.1  #
#  PAR1_Y8       6959       |    38106       #
#  PAR2_Y6       13971      |    30914       #
#  PAR3_y8       21188      |    22864       #
##############################################

PAR_data %>% filter(PAR > 10) %>%
  group_by(Group) %>% summarise(count = n()) %>%
  arrange(desc(count))

##### A tibble : 4 x 2#######    ######### 총 데이터 수 ###########
## Group        10 < X      #    #    Group         Data          #
#  PAR1_Y8       31         #    #   PAR1_Y8        45096         #
#  PAR2_Y6       211        #    #   PAR2_Y6        45096         #
#  PAR3_y8       1044       #    #   PAR3_y8        45096         #
#############################    ##################################