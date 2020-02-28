## source("~/Projects/coloc-cond-mask/load-hapdata.R")
setwd("~/D") 
library(randomFunctions)
library(magrittr)
library(data.table)
library(progress)
library(viridis)
library(ggplot2)
## devtools::load_all("~/RP/coloc")

dkeys <- "/home/cew54/share/Projects/twas/keys"
ddata <- sub("keys","sims",dkeys)  %>% paste0(.,"/results")
keyfiles <- list.files(dkeys,pattern="^simv6")
datafiles <- list.files(ddata,pattern="^simv6")
files <- intersect(keyfiles,datafiles)
message("keyfiles found: ",length(keyfiles))
message("matching output files found: ",length(files))
res <- vector("list",length(files))
pb <- progress_bar$new(total = length(files))
for(i in seq_along(files)) {
    pb$tick()
    f <- files[i]
    key <- readRDS(file.path(dkeys,f))
    data <- readRDS(file.path(ddata,f))  %>% as.data.table()
    data[,c("t","e1","e4","e2"):=key[c("t","e1","e4","e2")]]
    res[[i]] <- copy(data)
}
res %<>% rbindlist()
saveRDS(res, file = '/home/cew54/share/Projects/twas/sim_results.rds')


