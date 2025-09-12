###################
# 패턴별 코사인 유사도 결과 파일에 패턴 열 추가
###################

library(readr)    # read_csv / write_csv
library(dplyr)    # mutate
library(stringr)  # str_extract
library(purrr)    # walk

getwd()

setwd("C:/Users/USER/Desktop/Analysis/p-value/input/cosine_0.9_target/")

Sys.glob("*")


# 1) 대상 파일 목록 가져오기 ---------------------------------------------
files <- list.files(
  pattern = "^cosine_0\\.9_targets_\\(.*\\)\\.csv$",   # 괄호·점 이스케이프
  full.names = TRUE
)

# 2) 파일마다 pattern 열 추가 → 덮어쓰기 ----------------------------------
walk(files, function(f) {
  pat <- str_extract(basename(f), "\\(.*\\)")      # "(a,b,c)" 부분만 추출
  
  read_csv(f, show_col_types = FALSE) %>% 
    mutate(pattern = pat, .after = Name) %>%       # Name 뒤에 pattern 열 삽입
    write_csv(f)                                   # 같은 이름으로 저장
})

message(length(files), " files updated ✔")

#####################
# 코사인 유사도 결과와 logFC2, p-value 결과와 오버랩되는 이온 추출
#####################

library(readr)
library(dplyr)
library(stringr)
library(purrr)

setwd("../")
# 0) logFC·p-value 파일 읽기 -------------------------------------------------
df_log <- read_csv("ion_pattern_p-val_and_logFC_fit3_limma.csv", show_col_types = FALSE)
# df_log: Ion, pattern, logFC_*, pval_*, adjp_* …

# 1) 대상 패턴 파일 목록 ------------------------------------------------------
pat_files <- list.files(
  pattern = "^cosine_0\\.9_targets_\\(.*\\)\\.csv$",
  full.names = TRUE
)

# 2) 패턴별 루프 --------------------------------------------------------------
walk(pat_files, function(f) {
  
  # 2-1) 파일명에서 패턴 문자열 추출 例 "(1,0,1)"
  pat <- str_extract(basename(f), "\\(.*\\)")
  
  # 2-2) 패턴 파일 읽기
  df_pat <- read_csv(f, show_col_types = FALSE)      # pattern, Name, …  또는 빈 tibble
  
  # 2-3) 내용이 없으면 건너뜀
  if (nrow(df_pat) == 0) {
    message("skip (empty): ", basename(f))
    return(invisible(NULL))
  }
  
  # 2-4) logFC 데이터에서 같은 패턴만 선택
  df_log_pat <- df_log %>% filter(pattern == pat)
  if (nrow(df_log_pat) == 0) {
    message("skip (logFC 쪽에 해당 패턴 없음): ", pat)
    return(invisible(NULL))
  }
  
  # 2-5) Name = Ion & pattern 동시 일치 행만 추출
  overlap <- inner_join(
    df_pat,
    df_log_pat,
    by = c("Name" = "Ion", "pattern" = "pattern")
  )
  if (nrow(overlap) == 0) {
    message("skip (교집합 0): ", pat)
    return(invisible(NULL))
  }
  
  # 2-6) 결과 저장
  out_name <- paste0("overlap_", gsub("[()]", "", pat), ".csv")  # 1,0,1 등
  write_csv(overlap, out_name)
  message("saved: ", out_name, " (", nrow(overlap), " rows)")
})

# skip (empty): … -> 해당 cosine_0.9_targets_(a,b,c).csv 파일이 빈 파일이라서(행이 0개) 비교를 건너뜀.
# skip (logFC 쪽에 해당 패턴 없음): (a,b,c) -> ion_pattern_p-val_and_logFC.csv 안에 그 패턴값이 들어 있는 행 자체가 없음.
# skip (교집합 0): (a,b,c) -> 두 파일 모두에 (a,b,c) 패턴은 있지만, 이름(Name = Ion)이 겹치는 이온이 없어서 교집합이 0개.
# saved: overlap_a,b,c.csv (n rows) -> 두 파일에서 패턴과 이온 이름이 동시에 일치하는 행이 n개 나와서 overlap_a,b,c.csv로 저장 완료


################ 
#오버랩을 한 파일로 합치기
################ 
library(readr)
library(dplyr)
library(purrr)

# 1) overlap_* 파일 목록 수집 -------------------------------------------
overlap_files <- list.files(
  pattern = "^overlap_.*\\.csv$",    # overlap_1,0,1.csv 등
  full.names = TRUE
)

if (length(overlap_files) == 0) {
  stop("overlap_*.csv 파일이 없습니다. 먼저 교집합 파일을 생성하세요.")
}

# 2) 파일들을 차례로 읽어 행 방향으로 이어붙이기 --------------------------
df_all_overlap <- map_dfr(overlap_files, read_csv, show_col_types = FALSE)

# 3) (선택) 중복 이름·컬럼 순서 확인 후 저장 -----------------------------
# 필요에 따라 select()로 열 순서를 정렬하거나 distinct()로 행 중복 제거 가능
# df_all_overlap <- df_all_overlap %>% distinct()

write_csv(df_all_overlap, "overlap_all_patterns.csv")
message("✔  ", nrow(df_all_overlap), " rows가 'overlap_all_patterns.csv'에 저장되었습니다.")
