## source("~/Projects/coloc-cond-mask/load-hapdata.R")
library(randomFunctions)
library(magrittr)
library(data.table)
library(progress)
library(viridis)
library(ggplot2)
## devtools::load_all("~/RP/coloc")

##' make character vector from character string
ss <- function(x) strsplit(x,"")[[1]]

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
    data[,c("t","e1","e2","rho"):=key[c("t","e1","e2","rho")]]
    ## fix bug on storing e4
    data[,e4:=ifelse(key$e2=="-","-",
              ifelse(key$e2=="t",key$t,
              ifelse(key$e2=="same",key$e1,
                     setdiff(c("A","B","C","D"), c(ss(key$t),ss(key$e1)))[1])))]
    data[,LD.t.e1:=max(key$LD[ ss(key$t),ss(key$e1) ]^2)]
    data[,LD.t.e4:=max(key$LD[ ss(key$t),ss(data$e4[1]) ]^2)]
    data[,LD.e1.e4:=max(key$LD[ ss(key$e1),ss(data$e4[1]) ]^2)]
    ## data[,beta.t:=max(abs(key$beta[ key$t ]))]
    ## data[,beta.e1:=max(abs(key$beta[ key$e1 ]))]
    res[[i]] <- copy(data[trait=="E1"])
    ## res[[i]] <- copy(data)
}
res %<>% rbindlist(.,fill=TRUE)
saveRDS(res, file = '/home/cew54/share/Projects/twas/sim_results.rds')


