getwd()

setwd("./input/")

Sys.glob("*")

install.packages("dplyr")
install.packages("ggplot2")

library(dplyr)

data = data_frame(read.csv("Transition_Results_250527.csv"))

str(data)

m1_data <- data %>% filter(Method == 1)
m2_data <- data %>% filter(Method == 2)
m3_data <- data %>% filter(Method == 3)

m1_data <- arrange(m1_data, Replicate)
m2_data <- arrange(m2_data, Replicate)
m3_data <- arrange(m3_data, Replicate)

m1_IS <- m1_data %>% filter(grepl("^Peptide", Protein))
m1_patient <- m1_data %>% filter(!grepl("^Peptide", Protein))

m2_IS <- m2_data %>% filter(grepl("^Peptide", Protein))
m2_patient <- m2_data %>% filter(!grepl("^Peptide", Protein))

m3_IS <- m3_data %>% filter(grepl("^Peptide", Protein))
m3_patient <- m3_data %>% filter(!grepl("^Peptide", Protein))

prefix_nums <- c(20, 30, 40, 50)
suffix_nums <- c(1:6)  

data_list <- list(m1_data = m1_data, m2_data = m2_data, m3_data = m3_data)

all_data <- list()

for (name in names(data_list)) {
  df <- data_list[[name]]
  split_list <- list()
  
  for (i in prefix_nums) {
    for (j in suffix_nums) {
      key <- paste0(i, "s-", j)
      varname <- paste0("Rep_", i, "s_", j)
      
      subset_df <- df %>% filter(Replicate == key)
      
      if (nrow(subset_df) > 0) {
        split_list[[varname]] <- subset_df
      }
    }
  }
  
  all_data[[name]] <- split_list
}  

for (top_name in names(all_data)) {
  sublist <- all_data[[top_name]]
  
  for (sub_name in names(sublist)) {
    var_name <- paste0(top_name, "_", sub_name)
    
    assign(var_name, sublist[[sub_name]])
  }
}


extract_peptide_area <- function(df, df_name, peptide_id = "Peptide-03", n = 6) {
  target_rows <- df %>%
    filter(Protein == peptide_id) %>%
    head(n)
  
  for (i in 1:nrow(target_rows)) {
    varname <- paste0(df_name, "_stand", i)
    assign(varname, target_rows$Area[i], envir = .GlobalEnv)
  }
}

extract_peptide_area(m1_data_Rep_20s_1, "m1_20s_1")
extract_peptide_area(m1_data_Rep_20s_2, "m1_20s_2")
extract_peptide_area(m1_data_Rep_20s_3, "m1_20s_3")
extract_peptide_area(m1_data_Rep_20s_4, "m1_20s_4")
extract_peptide_area(m1_data_Rep_20s_5, "m1_20s_5")
extract_peptide_area(m1_data_Rep_20s_6, "m1_20s_6")

extract_peptide_area(m1_data_Rep_30s_1, "m1_30s_1")
extract_peptide_area(m1_data_Rep_30s_2, "m1_30s_2")
extract_peptide_area(m1_data_Rep_30s_3, "m1_30s_3")
extract_peptide_area(m1_data_Rep_30s_4, "m1_30s_4")
extract_peptide_area(m1_data_Rep_30s_5, "m1_30s_5")
extract_peptide_area(m1_data_Rep_30s_6, "m1_30s_6")

extract_peptide_area(m1_data_Rep_40s_1, "m1_40s_1")
extract_peptide_area(m1_data_Rep_40s_2, "m1_40s_2")
extract_peptide_area(m1_data_Rep_40s_3, "m1_40s_3")
extract_peptide_area(m1_data_Rep_40s_4, "m1_40s_4")
extract_peptide_area(m1_data_Rep_40s_5, "m1_40s_5")
extract_peptide_area(m1_data_Rep_40s_6, "m1_40s_6")

extract_peptide_area(m1_data_Rep_50s_1, "m1_50s_1")
extract_peptide_area(m1_data_Rep_50s_2, "m1_50s_2")
extract_peptide_area(m1_data_Rep_50s_3, "m1_50s_3")
extract_peptide_area(m1_data_Rep_50s_4, "m1_50s_4")
extract_peptide_area(m1_data_Rep_50s_5, "m1_50s_5")
extract_peptide_area(m1_data_Rep_50s_6, "m1_50s_6")



extract_peptide_area(m2_data_Rep_20s_1, "m2_20s_1")
extract_peptide_area(m2_data_Rep_20s_2, "m2_20s_2")
extract_peptide_area(m2_data_Rep_20s_3, "m2_20s_3")
extract_peptide_area(m2_data_Rep_20s_4, "m2_20s_4")
extract_peptide_area(m2_data_Rep_20s_5, "m2_20s_5")
extract_peptide_area(m2_data_Rep_20s_6, "m2_20s_6")

extract_peptide_area(m2_data_Rep_30s_1, "m2_30s_1")
extract_peptide_area(m2_data_Rep_30s_2, "m2_30s_2")
extract_peptide_area(m2_data_Rep_30s_3, "m2_30s_3")
extract_peptide_area(m2_data_Rep_30s_4, "m2_30s_4")
extract_peptide_area(m2_data_Rep_30s_5, "m2_30s_5")
extract_peptide_area(m2_data_Rep_30s_6, "m2_30s_6")

extract_peptide_area(m2_data_Rep_40s_1, "m2_40s_1")
extract_peptide_area(m2_data_Rep_40s_2, "m2_40s_2")
extract_peptide_area(m2_data_Rep_40s_3, "m2_40s_3")
extract_peptide_area(m2_data_Rep_40s_4, "m2_40s_4")
extract_peptide_area(m2_data_Rep_40s_5, "m2_40s_5")
extract_peptide_area(m2_data_Rep_40s_6, "m2_40s_6")

extract_peptide_area(m2_data_Rep_50s_1, "m2_50s_1")
extract_peptide_area(m2_data_Rep_50s_2, "m2_50s_2")
extract_peptide_area(m2_data_Rep_50s_3, "m2_50s_3")
extract_peptide_area(m2_data_Rep_50s_4, "m2_50s_4")
extract_peptide_area(m2_data_Rep_50s_5, "m2_50s_5")
extract_peptide_area(m2_data_Rep_50s_6, "m2_50s_6")


 
extract_peptide_area(m3_data_Rep_20s_1, "m3_20s_1")
extract_peptide_area(m3_data_Rep_20s_2, "m3_20s_2")
extract_peptide_area(m3_data_Rep_20s_3, "m3_20s_3")
extract_peptide_area(m3_data_Rep_20s_4, "m3_20s_4")
extract_peptide_area(m3_data_Rep_20s_5, "m3_20s_5")
extract_peptide_area(m3_data_Rep_20s_6, "m3_20s_6")

extract_peptide_area(m3_data_Rep_30s_1, "m3_30s_1")
extract_peptide_area(m3_data_Rep_30s_2, "m3_30s_2")
extract_peptide_area(m3_data_Rep_30s_3, "m3_30s_3")
extract_peptide_area(m3_data_Rep_30s_4, "m3_30s_4")
extract_peptide_area(m3_data_Rep_30s_5, "m3_30s_5")
extract_peptide_area(m3_data_Rep_30s_6, "m3_30s_6")

extract_peptide_area(m3_data_Rep_40s_1, "m3_40s_1")
extract_peptide_area(m3_data_Rep_40s_2, "m3_40s_2")
extract_peptide_area(m3_data_Rep_40s_3, "m3_40s_3")
extract_peptide_area(m3_data_Rep_40s_4, "m3_40s_4")
extract_peptide_area(m3_data_Rep_40s_5, "m3_40s_5")
extract_peptide_area(m3_data_Rep_40s_6, "m3_40s_6")

extract_peptide_area(m3_data_Rep_50s_1, "m3_50s_1")
extract_peptide_area(m3_data_Rep_50s_2, "m3_50s_2")
extract_peptide_area(m3_data_Rep_50s_3, "m3_50s_3")
extract_peptide_area(m3_data_Rep_50s_4, "m3_50s_4")
extract_peptide_area(m3_data_Rep_50s_5, "m3_50s_5")
extract_peptide_area(m3_data_Rep_50s_6, "m3_50s_6")

# stand1 : YVYVADVAAK -> y8 (top transition)
# stand2 : YVYVADVAAK -> y6 (top transition)
# stand3 : YVYVADVAAK -> y8 (2nd transition)
# stand4 : YVYVADVAAK -> y6 (2nd transition)
# stand5 : YVYVADVAAK -> y8 (3rd transition)
# stand6 : YVYVADVAAK -> y6 (3rd transition)

save.image("./Aging_Clock_data_preprocessing_250605.RData")
load("./Aging_Clock_data_preprocessing_250605.RData")


generate_PAR_columns <- function(df_prefix = "m1_data_Rep", stand_prefix = "m1", 
                                 prefix_nums = c(20, 30, 40, 50), suffix_nums = 1:6) {
  
  for (pre in prefix_nums) {
    for (suf in suffix_nums) {
      df_name <- paste0(df_prefix, "_", pre, "s_", suf)
      stand_base <- paste0(stand_prefix, "_", pre, "s_", suf)
      
      if (exists(df_name, envir = .GlobalEnv)) {
        df <- get(df_name, envir = .GlobalEnv)
        
        for (i in 1:6) {
          stand_var <- paste0(stand_base, "_stand", i)
          if (exists(stand_var, envir = .GlobalEnv)) {
            stand_val <- get(stand_var, envir = .GlobalEnv)
            
            if (!is.na(stand_val) && stand_val != 0) {
              df[[paste0("PAR", i)]] <- df$Area / stand_val
            } else {
              df[[paste0("PAR", i)]] <- NA
            }
          }
        }
        
        assign(df_name, df, envir = .GlobalEnv)
      }
    }
  }
}

generate_PAR_columns()

generate_PAR_columns(df_prefix = "m2_data_Rep", stand_prefix = "m2")   # Method 2

generate_PAR_columns(df_prefix = "m3_data_Rep", stand_prefix = "m3")   # Method 3

combined_data <- function(data_prefix = "m1_data_Rep",
                             prefix_nums = c(20, 30, 40, 50),
                             suffix_nums = 1:6) {
  
  df_list <- list()
  
  for (pre in prefix_nums) {
    for (suf in suffix_nums) {
      df_name <- paste0(data_prefix, "_", pre, "s_", suf)
                      #m1_data_Rep_20s_01
      
      if (exists(df_name, envir = .GlobalEnv)) {
        df <- get(df_name, envir = .GlobalEnv)
        df_list[[df_name]] <- df
      }
    }
  }
  
  combined_df <- do.call(rbind, df_list)
  return(combined_df)
}

m1_df <- combined_data(data_prefix="m1_data_Rep")
m2_df <- combined_data(data_prefix="m2_data_Rep")
m3_df <- combined_data(data_prefix="m3_data_Rep")

write.csv(m1_df, file="../output/m1_result_PAR.csv", row.names=F)
write.csv(m2_df, file="../output/m2_result_PAR.csv", row.names=F)
write.csv(m3_df, file="../output/m3_result_PAR.csv", row.names=F)

final_df <- rbind(m1_df, m2_df)
final_df <- rbind(final_df, m3_df)

write.csv(final_df, file="../output/Transition_PAR_All (IS included)_250527.csv", row.names=F)

final_IS <- final_df %>% filter(grepl("^Peptide", Protein))
final_Sample <- final_df %>% filter(!grepl("^Peptide", Protein))  # Keep only the sample peptides (exclude internal standards)

write.csv(final_Sample, file="../output/Transition_PAR_Results_250605.csv", row.names=F)

collect_stand_values <- function(stand_prefix = "m1", 
                                 prefix_nums = c(20, 30, 40, 50), 
                                 suffix_nums = 1:6) {
  
  rows_list <- list()
  
  for (pre in prefix_nums) {
    for (suf in suffix_nums) {
      id <- paste0(pre, "s-", suf)
      base_name <- paste0(stand_prefix, "_", pre, "s_", suf)
      
      row_vals <- numeric(3)
      complete <- TRUE
      
      for (i in 1:6) {
        var_name <- paste0(base_name, "_stand", i)
        if (exists(var_name, envir = .GlobalEnv)) {
          row_vals[i] <- get(var_name, envir = .GlobalEnv)
        } else {
          row_vals[i] <- NA
          complete <- FALSE
        }
      }
      
      rows_list[[id]] <- row_vals
    }
  }
  

  stand_df <- do.call(rbind, rows_list) %>% as.data.frame()
  colnames(stand_df) <- c("stand1", "stand2", "stand3", "stand4", "stand5", "stand6")
  stand_df$Replicate <- rownames(stand_df)
  rownames(stand_df) <- NULL
  
  stand_df <- stand_df[, c("Replicate", "stand1", "stand2", "stand3", "stand4", "stand5", "stand6")]
  
  return(stand_df)
}

m1_stand_table <- collect_stand_values(stand_prefix="m1")
m2_stand_table <- collect_stand_values(stand_prefix="m2")
m3_stand_table <- collect_stand_values(stand_prefix="m3")

stand_data <- rbind(m1_stand_table, m2_stand_table)
stand_data <- rbind(stand_data, m3_stand_table)

m_list <- c(rep(1, 24), rep(2, 24), rep(3, 24))
stand_data$method <- m_list

stand_data <- stand_data %>% relocate(method, .before = stand1)

write.csv(stand_data, file="../output/Transition_IS_area_250605.csv", row.names=F)


### Figure 1B bar plot. Fraction of PAR values within the target range (0.1 ≤ PAR ≤ 10)
library(ggplot2)
library(dplyr)
library(scales)

# bar plot (final)
ggplot(band_rate_all, aes(x = stand, y = band_rate, fill = stand)) +
  geom_col(width = 0.6, color = "gray40") +  # bar width / border
  geom_text(
    aes(label = percent(band_rate, accuracy = 0.1)),
    vjust = -0.4, size = 4, color = "black"
  ) +
  scale_fill_manual(
    name = "Internal Standard",
    values = c(
      "stand1" = "#F4A7B9",  # light apricot pink (PAR1_Y8)
      "stand2" = "#D6C26E",  # mustard beige (PAR2_Y6)
      "stand3" = "#8FD19E",  # mint green (PAR3_Y8)
      "stand4" = "#7DCFE0",  # sky blue (PAR4_Y6)
      "stand5" = "#8DA9E0",  # blue-gray (PAR5_Y8)
      "stand6" = "#D79ACD"   # lilac pink (PAR6_Y6)
    )
  ) +
  labs(
    x = "Internal Standard (stand)",
    y = "Fraction within 0.1–10 (log10 −1~1)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",               # place color legend on the right
    legend.title = element_text(face = "bold"),
    axis.text = element_text(color = "gray25"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )
# band_rate(stand) = N_in(stand) / N_total(stand) -> rate = (count of valid observations where logPAR is not NA) / (count of rows satisfying -1 <= logPAR <= 1)
# How stably each internal standard (stand) normalized the PAR -> the key metric used as the basis for deciding which IS to select
