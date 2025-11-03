library(readr)
library(missMDA)
library(dplyr)

sep_df <- read_csv("data/df_analysis.csv",
             col_types = cols(
               Subject_ID = col_character(),
               group = col_character(),
               time = col_double(),
               errors = col_double(),
               front = col_double(),
               all = col_double(),
               speed = col_double(),
               nboli = col_double(),
               polig = col_double(),
               micro = col_double(),
               mog = col_double(),
               ftl = col_double(),
               pmog = col_double(),
               pftl = col_double()
             ))

group <- sep_df$group
subject <- sep_df$Subject_ID


sep_df <- sep_df %>% select(-group)
sep_df <- sep_df %>% select(-Subject_ID)
sep_df <- sep_df %>% select(-nboli)


nb <- estim_ncpPCA(sep_df, ncp.min=0, ncp.max=5)
res_impute <- imputePCA(sep_df, ncp = nb$ncp)
sep_df <- as.data.frame(res_impute$completeObs)

col_clinique <- c('time', 'errors', 'front', 'all', 'speed')
col_biologique <- c('polig', 'micro', 'mog', 'ftl', 'pmog', 'pftl')

