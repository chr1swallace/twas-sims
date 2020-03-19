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

BATCH.SIZE <- 20 # how many files to process per go
taskid <- 0
library(randomFunctions)
args <- getArgs()
print(args)
if("BATCH.SIZE" %in% names(args))
  BATCH.SIZE <- as.numeric(args$BATCH.SIZE)
if("taskid" %in% names(args))
  taskid=as.numeric(args$taskid)
  

setwd("/home/cew54/share/Projects/twas/sims")
do.gwas <- function(y, geno, stratum = NULL, cc = FALSE) {
  res <- lapply(1 : ncol(geno), FUN = function(j) {
    ## message(colnames(geno)[j])
    m <- lm(y ~ geno[, j])
   summ <- summary(m) 
    data.frame(snp = colnames(geno)[j],
               beta = summ$coefficients[2, "Estimate"],
               var.beta = summ$coefficients[2, "Std. Error"]^2,
               p.value = summ$coefficients[2, "Pr(>|t|)"])
  })
  do.call(rbind, res)
}

FF <- list.files(pattern = "simv6")
## FF <- FF[!(FF %in% err$file)]
files.done <- list.files("results")
files.todo <- setdiff(FF, files.done)[ (taskid * BATCH.SIZE + 1):((taskid+1)*BATCH.SIZE) ]   
message("running ",length(files.todo)," files.")

for(m in sample(files.todo)) {
  if(file.exists(file.path("results", m)))
    next
  message("\n",m)
  x <- readRDS(m)
  y <- x$y
  z <- x$z; colnames(z) <- paste0("E", 1 : ncol(z))
  Gy <- as(x[[2]], "numeric")
  Gz <- as(x[[4]], "numeric")
  pocket <- names(x$pocket)

  ##QC for monomorphic SNPs
  ## library(microbenchmark)
  ## microbenchmark(ind.z <- apply(Gz, 2, FUN = function(u) length(unique(u)) > 1),
  ##                ind.z <- apply(Gz, 2, sd) == 0)
  ind.y <- apply(Gy, 2, sd) != 0
  ind.z <- apply(Gz, 2, sd) != 0
  ind <- ind.z & ind.y
  if(any(!ind)) {
    Gz <- Gz[, ind]
    Gy <- Gy[, ind]
    if(length(pocket))
      pocket <- intersect(pocket,colnames(Gz)) #make sure get rid of pocket SNPs that were removed through QC
  }
    
  ## z.mtl <- scale(reshape2::melt(z)$value)
  ## z.mtl <- as.vector(scale(z))
  z.mtl2 <- scale(as.vector(z))
# y.mtl <- scale(c(y, y, y, y, y))
# Gz.mtl <- do.call(rbind, list(Gz, Gz, Gz, Gz, Gz))
# Gy.mtl <- do.call(rbind, list(Gy, Gy, Gy, Gy, Gy))
  y.mtl <- scale(rep(y, ncol(z)))
  Gz.mtl <- Gy.mtl <- vector('list', ncol(z))
  for(i in 1 : ncol(z)) {
    Gz.mtl[[i]] <- Gz
    Gy.mtl[[i]] <- Gy
  }
  Gz.mtl <- do.call(rbind, Gz.mtl)
  Gy.mtl <- do.call(rbind, Gy.mtl)
  id <- rep(colnames(z), each = nrow(z))
  z <- scale(z)
  y <- scale(y)
#---------------------------------------

##GWAS
# snp.info <- strsplit(colnames(Gz), split = "\\.")
# snp.info <- data.frame(snp = sapply(snp.info, "[[", 1), pos = as.numeric(sapply(snp.info, "[[", 2)))
z.gwas <- apply(z, 2, FUN = function(w) do.gwas(w, Gz))
y.gwas <- do.gwas(y, Gy)

# par(mfrow = c(2, 3))
# for(i in 1 : 5) plot(snp.info$pos, -log10(z.gwas[[i]]$p.value), main = colnames(z)[i])
# plot(snp.info$pos, -log10(y.gwas$p.value), main = "GWAS trait")

##LASSO
l.preds <- lapply(1 : ncol(z), FUN = function(i) {
  message(sprintf("LASSO: starting %s...", colnames(z)[i]))

  mod <- cv.glmnet(Gz, z[, i], family = "gaussian", alpha = 1, standardize = FALSE, parallel = FALSE)
  pred <- predict(mod, Gy)
  message(sprintf("LASSO: %s finished!", colnames(z)[i]))

  pred
}) %>% do.call(cbind, .)

## system.time(l.preds <- lapply(1 : ncol(z), FUN = function(i) {
##   message(sprintf("LASSO: starting %s...", colnames(z)[i]))

##   mod <- cv.glmnet(Gz, z[, i], family = "gaussian", alpha = 1, standardize = FALSE, parallel = TRUE)
##   pred <- predict(mod, Gy)
##   message(sprintf("LASSO: %s finished!", colnames(z)[i]))

##   pred
## }) %>% do.call(cbind, .))

##RF
rf.preds <- lapply(1 : ncol(z), FUN = function(i) {
  message(sprintf("RF: starting %s...", colnames(z)[i]))

  mod <- ranger(y ~ ., data.frame(Gz, y = z[, i]), num.trees = 500, mtry = floor(ncol(Gz)/3), importance = "permutation")
  pred <- predict(mod, Gy)$prediction
  message(sprintf("RF: %s finished!", colnames(z)[i]))

  pred
}) %>% do.call(cbind, .)

##JOINT LASSO
lambda <- exp(seq(-20, 0, length = 20))

corr <- cor(z, use = 'pairwise.complete.obs')
G <- abs(corr)
diag(G) <- 0
G <- G/max(G)

## fu.mod <- L2fuse.cv(Gz.mtl, z.mtl, id, G, lambda)
## fu.pred <- L2fuse.pred(fu.mod$opt.beta, Gy.mtl, id)[, 1]
## fu.preds <- matrix(fu.pred,ncol=ncol(z)) #split(fu.pred, ceiling(seq_along(fu.pred)/nrow(z))) %>% do.call(cbind, .)


fu.mod2 <- L2fuse.cv(Gz.mtl, z.mtl2, id, G, lambda)
fu.pred2 <- L2fuse.pred(fu.mod2$opt.beta, Gy.mtl, id)[, 1]
fu.preds2 <- matrix(fu.pred2,ncol=ncol(z)) #split(fu.pred, ceiling(seq_along(fu.pred)/nrow(z))) %>% do.call(cbind, .)

#RF-MTL
  ## mtl.mod <- ranger(y ~ ., data = data.frame(Gz.mtl, id = id, y = z.mtl), num.trees = 500, mtry = floor(ncol(Gz.mtl)/3), importance = "permutation", always.split.variables = "id")
  ## mtl.pred <- predict(mtl.mod, data.frame(Gy.mtl, id))$predictions
  ## mtl.preds <- matrix(mtl.pred,ncol=ncol(z)) # split(mtl.pred, ceiling(seq_along(mtl.pred)/nrow(z))) %>% do.call(cbind, .)

  mtl.mod2 <- ranger(y ~ ., data = data.frame(Gz.mtl, id = id, y = z.mtl2), num.trees = 500, mtry = floor(ncol(Gz.mtl)/3), importance = "permutation", always.split.variables = "id")
  mtl.pred2 <- predict(mtl.mod2, data.frame(Gy.mtl, id))$predictions
  mtl.preds2 <- matrix(mtl.pred2,ncol=ncol(z)) # split(mtl.pred, ceiling(seq_along(mtl.pred)/nrow(z))) %>% do.call(cbind, .)

##TWAS-association
do.twas <- function(preds, y) {
  apply(preds, 2, FUN = function(w) {
    coefs <- summary(lm(y ~ w))$coefficients
    if(nrow(coefs) == 2) return(coefs[2, 4]) else return(NA)
  })
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



res <- data.frame(trait = colnames(z),
  rf.twas = do.twas(rf.preds, y),
  lasso.twas = do.twas(l.preds, y),
  mtl2.twas = do.twas(mtl.preds2, y),
  fuser2.twas = do.twas(fu.preds2, y),
  ## mtl.twas = do.twas(mtl.preds, y),
  ## fuser.twas = do.twas(fu.preds, y),
  coloc.pval = coloc.res,
  min.pval.expr = sapply(z.gwas, FUN = function(w) min(w$p.value)),
  min.snp.expr = sapply(z.gwas, FUN = function(w) w$snp[which.min(w$p.value)]),
  min.pval.gwas = min(y.gwas$p.value),
  min.snp.gwas = y.gwas$snp[which.min(y.gwas$p.value)])


  saveRDS(res, file = file.path("results/", m))
  message(sprintf("%s: done!", m))
}
