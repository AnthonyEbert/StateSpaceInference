
library(dplyr)
x <- data.frame(value = NA, weight = NA, seed = NA, type = NA)

for(i in 0:14){
  x <- x %>%
    dplyr::bind_rows(readRDS(paste0("theta_df_", i, ".RData"))) %>%
    dplyr::bind_rows(readRDS(paste0("theta_stan_", i, ".RData")))
}

saveRDS(x, file = "combined_df.rds")
