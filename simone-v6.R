setwd("~/D")
source("~/Projects/coloc-cond-mask/load-hapdata.R")
library(glmnet)
library(fuser)
library(ranger)
library(annotSnpStats)
library(coloc)
## library(doMC)
## registerDoMC(cores = 2)
source("fuser_hack_nfg.R")
source("Code_util.R")
library(magrittr)
library(parallel)
source("common.R")

## now we add 2 cvs
library(randomFunctions)
library(magrittr)
devtools::load_all("~/RP/coloc")
args <- getArgs(defaults=list(N=400,NSIM=2,NSNP=500,ld="highld",pop="eur",
e1="A", # A, AB, -
e2="diff", # same, weaker, diff, -
t="A",rho=0), # A, AB, BA, AC, C
                numeric=c("N","NSIM","NCV","rho"))
print(args)

## hdata <- readRDS(paste0("~/scratch/",args$ld,".RDS"))
GG <- readRDS("/home/cew54/share/Projects/twas/sample_geno.rds")

do.twas <- function(data) {
  x <- data$data
  y <- x$y
  z <- x$z; colnames(z) <- paste0("E", 1 : ncol(z))
  Gy <- as(x[[2]], "numeric")
  Gz <- as(x[[4]], "numeric")
  pocket <- names(x$pocket)

  id <- rep(colnames(z), each = nrow(z))
  z <- scale(z[,1],scale=FALSE,center=TRUE)
  y <- scale(y,center=TRUE,scale=FALSE)
#---------------------------------------

##GWAS
# snp.info <- strsplit(colnames(Gz), split = "\\.")
# snp.info <- data.frame(snp = sapply(snp.info, "[[", 1), pos = as.numeric(sapply(snp.info, "[[", 2)))
  z.gwas <- do.gwas(z, Gz)
  y.gwas <- do.gwas(y, Gy)

##LASSO

  mod <- cv.glmnet(Gz, z[, 1], family = "gaussian", alpha = 1, standardize = FALSE, parallel = FALSE)
  lasso <- predict(mod, Gy)
  message(sprintf("LASSO: finished!"))

##RF
  mod <- ranger(y ~ ., data.frame(Gz, y = z[, 1]), num.trees = 500, mtry = floor(ncol(Gz)/3), importance = "permutation")
  rf <- predict(mod, Gy)$prediction
  message(sprintf("RF: finished!"))

##TWAS-association
twas <- function(w, y) {
  m <- lm(y ~ w)
  coefs <- summary(m)$coefficients
  ret <- if(nrow(coefs) == 2) coefs[2, c(1,2,4)] else rep(NA,3)
  names(ret) <- c("est","se","p")
  ret
}

#------------------------
##coloc
if(length(pocket)>1) {
  pcs <- pcs.prepare(Gz[, pocket], Gy[, pocket])

  npcs <- min(which(pcs@vars > 0.8)[1], 6)
  if(npcs < 2) npcs <- 2
  thr <- pcs@vars[npcs - 1]
  m1 <- pcs.model(pcs, group = 1, Y = z[,1], family = "gaussian", threshold = thr)
  m2 <- pcs.model(pcs, group = 2, Y = y, family = "gaussian", threshold = thr)
  tmp <- coloc.test(m1, m2, bayes = FALSE)
  
  chisqr <- tmp@result['chisquare']
  coloc.res <- pchisq(chisqr, tmp@result['n'] - 1, lower.tail = FALSE)
} else {
  coloc.res <- as.numeric(NA)
}

res <- c(rf = twas(rf, y),
         lasso = twas(lasso, y),
         coloc.pval = coloc.res)  %>% as.list()  %>% as.data.frame()
}

for(i in 1:args$NSIM) {
  data <- NULL
  while(is.null(data))
    data=simone()
  twas <- do.twas(data)

  ## dir.create("~/share/Projects/twas/sims",recursive=TRUE)
  ## dir.create("~/scratch/Projects/twas/sims",recursive=TRUE)
  
  tmp <- tempfile(pattern="simv6","~/share/Projects/twas/simone",fileext=".rds")
  key <- sub("sims","keys",tmp)
  
  ## saveRDS(x$data, file=tmp, version=2)
  saveRDS(list(key=x$key,twas=twas), file=key, version=2)
}


