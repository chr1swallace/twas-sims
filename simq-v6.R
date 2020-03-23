source("~/Projects/coloc-cond-mask/load-hapdata.R")
source("~/D/common.R")

## now we add 2 cvs
library(randomFunctions)
library(magrittr)
devtools::load_all("~/RP/coloc")
args <- getArgs(defaults=list(N=400,NSIM=2,NSNP=500,ld="highld",pop="eur",
e1="A", # A, AB, -
e2="same", # same, weaker, diff, -
t="A",rho=0), # A, AB, BA, AC, C
                numeric=c("N","NSIM","NCV","rho"))
print(args)

## hdata <- readRDS(paste0("~/scratch/",args$ld,".RDS"))
GG <- readRDS("/home/cew54/share/Projects/twas/sample_geno.rds")
for(i in 1:args$NSIM) {
  x <- NULL
  while(is.null(x))
    x=simone()
  ## dir.create("~/share/Projects/twas/sims",recursive=TRUE)
  ## dir.create("~/scratch/Projects/twas/sims",recursive=TRUE)
  
  tmp <- tempfile(pattern="simv6","~/share/Projects/twas/sims",fileext=".rds")
  key <- sub("sims","keys",tmp)
  
  saveRDS(x$data, file=tmp, version=2)
  saveRDS(x$key, file=key, version=2)
}


