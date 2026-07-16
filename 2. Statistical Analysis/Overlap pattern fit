library(readr)
library(dplyr)
library(stringr)
library(purrr)

getwd()

setwd("C:/Users/USER/Desktop/Analysis/p-value/input/cosine_0.9_target/")

Sys.glob("*")

files <- list.files(
  pattern = "^cosine_0\\.9_targets_\\(.*\\)\\.csv$",
  full.names = TRUE
)


walk(files, function(f) {
  pat <- str_extract(basename(f), "\\(.*\\)")
  
  read_csv(f, show_col_types = FALSE) %>% 
    mutate(pattern = pat, .after = Name) %>% 
    write_csv(f)
})

message(length(files), " files updated ✔")



library(readr)
library(dplyr)
library(stringr)
library(purrr)

setwd("../")

df_log <- read_csv("ion_pattern_p-val_and_logFC_fit3_limma.csv", show_col_types = FALSE)

pat_files <- list.files(
  pattern = "^cosine_0\\.9_targets_\\(.*\\)\\.csv$",
  full.names = TRUE
)


walk(pat_files, function(f) {
  
  pat <- str_extract(basename(f), "\\(.*\\)")
  
  df_pat <- read_csv(f, show_col_types = FALSE)
  
  if (nrow(df_pat) == 0) {
    message("skip (empty): ", basename(f))
    return(invisible(NULL))
  }
  
  df_log_pat <- df_log %>% filter(pattern == pat)
  if (nrow(df_log_pat) == 0) {
    message("skip (logFC 쪽에 해당 패턴 없음): ", pat)
    return(invisible(NULL))
  }
  
  overlap <- inner_join(
    df_pat,
    df_log_pat,
    by = c("Name" = "Ion", "pattern" = "pattern")
  )
  if (nrow(overlap) == 0) {
    message("skip (교집합 0): ", pat)
    return(invisible(NULL))
  }
  
  out_name <- paste0("overlap_", gsub("[()]", "", pat), ".csv")  # 1,0,1 등
  write_csv(overlap, out_name)
  message("saved: ", out_name, " (", nrow(overlap), " rows)")
})



library(readr)
library(dplyr)
library(purrr)


overlap_files <- list.files(
  pattern = "^overlap_.*\\.csv$",    # overlap_1,0,1.csv 등
  full.names = TRUE
)

if (length(overlap_files) == 0) {
  stop("overlap_*.csv 파일이 없습니다. 먼저 교집합 파일을 생성하세요.")
}

df_all_overlap <- map_dfr(overlap_files, read_csv, show_col_types = FALSE)


write_csv(df_all_overlap, "overlap_all_patterns.csv")
message("✔  ", nrow(df_all_overlap), " rows가 'overlap_all_patterns.csv'에 저장되었습니다.")
