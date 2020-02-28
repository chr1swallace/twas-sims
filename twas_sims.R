library(glmnet)
library(fuser)
library(ranger)
library(annotSnpStats)
library(coloc)
library(doMC)
registerDoMC(cores = 10)

setwd("/home/cew54/share/Projects/twas/sims")
source("/rds/user/ng414/hpc-work/new_twas/fuser_hack_nfg.R")
source("/rds/user/ng414/hpc-work/new_twas/fuser_stuff/Code_util.R")

if(!file.exists("error_log.rds")) {
	err <- data.frame(file = character(), error = character())
	saveRDS(err, file = "error_log.rds")
}
#----------------->>>>>FUNCTIONS

do.gwas <- function(y, geno, stratum = NULL, cc = FALSE) {
	if(cc) {

			res <- lapply(1 : ncol(geno), FUN = function(j) {
				message(colnames(geno)[j])
				if(!is.null(stratum)) {
					m <- glm(y ~ geno[, j] + stratum, family = 'binomial')
				} else {
						m <- glm(y ~ geno[, j], family = 'binomial')
				}
				data.frame(snp = colnames(geno)[j],
					beta = summary(m)$coefficients[2, "Estimate"],
					var.beta = summary(m)$coefficients[2, "Std. Error"]^2,
					p.value = summary(m)$coefficients[2, "Pr(>|z|)"])
			})

	} else {
		res <- lapply(1 : ncol(geno), FUN = function(j) {
			message(colnames(geno)[j])
			m <- lm(y ~ geno[, j])

			data.frame(snp = colnames(geno)[j],
				beta = summary(m)$coefficients[2, "Estimate"],
				var.beta = summary(m)$coefficients[2, "Std. Error"]^2,
				p.value = summary(m)$coefficients[2, "Pr(>|t|)"])
		})
	}

	do.call(rbind, res)
}
#----------------->>>>>FUNCTIONS
err <- readRDS("error_log.rds")

FF <- list.files(pattern = "simv5")
FF <- FF[!(FF %in% err$file)]

N <- length(list.files("results/", pattern = "simv5"))
if(N == length(FF)) stop("All simulations done!")

m <- sample(FF, 1)
while(file.exists(file.path("results", m))) {
	message(sprintf("%s already done, moving on...", m))
  m <- sample(FF, 1)
}
message(sprintf("doing: %s", m))

x <- readRDS(m)
y <- x$y
z <- x$z; colnames(z) <- paste0("E", 1 : ncol(z))
Gy <- as(x[[2]], "numeric")
Gz <- as(x[[4]], "numeric")
pocket <- names(x$pocket)
if(length(pocket) < 1) {
	err <- readRDS("error_log.rds")
	err <- rbind(err, data.frame(file = m, error = sprintf("only %s SNP in pocket", length(pocket))))
	err <- err[!duplicated(err$file), ]
	saveRDS(err, file = "error_log.rds")
}

#QC for monomorphic SNPs
ind.z <- apply(Gz, 2, FUN = function(u) length(unique(u)) > 1)
ind.y <- apply(Gy, 2, FUN = function(u) length(unique(u)) > 1)
ind <- ind.z & ind.y
Gz <- Gz[, ind]
Gy <- Gy[, ind]
pocket <- pocket[pocket %in% colnames(Gz)] #make sure get rid of pocket SNPs that were removed through QC

z.mtl <- scale(melt(z)$value)
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

  mod <- cv.glmnet(Gz, z[, i], family = "gaussian", alpha = 1, standardize = FALSE, parallel = TRUE)
  pred <- predict(mod, Gy)
  message(sprintf("LASSO: %s finished!", colnames(z)[i]))

  pred
}) %>% do.call(cbind, .)

##RF
rf.preds <- mclapply(1 : ncol(z), FUN = function(i) {
  message(sprintf("RF: starting %s...", colnames(z)[i]))

  mod <- ranger(y ~ ., data.frame(Gz, y = z[, i]), num.trees = 500, mtry = floor(ncol(Gz)/3), importance = "permutation")
  pred <- predict(mod, Gy)$prediction
  message(sprintf("RF: %s finished!", colnames(z)[i]))

  pred
}, mc.cores = 5) %>% do.call(cbind, .)

##JOINT LASSO
lambda <- exp(seq(-20, 0, length = 20))

corr <- cor(z, use = 'pairwise.complete.obs')
G <- abs(corr)
diag(G) <- 0
G <- G/max(G)

fu.mod <- L2fuse.cv(Gz.mtl, z.mtl, id, G, lambda)
fu.pred <- L2fuse.pred(fu.mod$opt.beta, Gy.mtl, id)[, 1]
fu.preds <- split(fu.pred, ceiling(seq_along(fu.pred)/nrow(z))) %>% do.call(cbind, .)

#RF-MTL
mtl.mod <- ranger(y ~ ., data = data.frame(Gz.mtl, id = id, y = z.mtl), num.trees = 500, mtry = floor(ncol(Gz.mtl)/3), importance = "permutation", always.split.variables = "id")
mtl.pred <- predict(mtl.mod, data.frame(Gy.mtl, id))$predictions
mtl.preds <- split(mtl.pred, ceiling(seq_along(mtl.pred)/nrow(z))) %>% do.call(cbind, .)


##TWAS-association
do.twas <- function(preds, y) {
  apply(preds, 2, FUN = function(w) {
    coefs <- summary(lm(y ~ w))$coefficients
    if(nrow(coefs) == 2) return(coefs[2, 4]) else return(NA)
  })
}

#------------------------
##coloc
pcs <- pcs.prepare(Gz[, pocket], Gy[, pocket])

coloc.res <- lapply(1 : ncol(z), FUN = function(i) {
  w <- z[, i]

  npcs <- min(which(pcs@vars > 0.8)[1], 6)
  if(npcs < 2) npcs <- 2
  thr <- pcs@vars[npcs - 1]

  m1 <- pcs.model(pcs, group = 1, Y = w, family = "gaussian", threshold = thr)
	m2 <- pcs.model(pcs, group = 2, Y = y, family = "gaussian", threshold = thr)
  tmp <- coloc.test(m1, m2, bayes = TRUE)
  chisqr <- tmp@result['chisquare']
  pchisq(chisqr, tmp@result['n'] - 1, lower.tail = FALSE)
}) %>% unlist




res <- data.frame(trait = colnames(z),
  rf.twas = do.twas(rf.preds, y),
  lasso.twas = do.twas(l.preds, y),
  mtl.twas = do.twas(mtl.preds, y),
  fuser.twas = do.twas(fu.preds, y),
  coloc.pval = coloc.res,
  min.pval.expr = sapply(z.gwas, FUN = function(w) min(w$p.value)),
  min.snp.expr = sapply(z.gwas, FUN = function(w) w$snp[which.min(w$p.value)]),
  min.pval.gwas = min(y.gwas$p.value),
  min.snp.gwas = y.gwas$snp[which.min(y.gwas$p.value)])


saveRDS(res, file = file.path("results/", m))
message(sprintf("%s: done!", m))

q('no')









##
#analysing results

FF <- list.files("results", full.names = TRUE)
info <- lapply(FF, FUN = function(a) {
  message(a)
  z <- readRDS(a)
  data.frame(z, id = gsub("results/|\\.rds", "", a))
}) %>% do.call(rbind, .)

info <- data.table(info)
ind <- rowSums(info[, 2 : 5] < 0.05, na.rm = TRUE) > 0
info.signif <- info[ind & coloc.pval > 0.05, ]

saveRDS(info.signif, file = "signif_ress.rds")
saveRDS(info, file = "results.rds")



"/rds/project/cew54/rds-cew54-wallace-share/Projects/twas/sims/signif_ress.rds"


























##

#
# coloc.res <- lapply(1 : ncol(z), FUN = function(i) {
#   message(i)
#   w <- z[, i]
#
#   ind.y <- y.gwas[y.gwas$p.value < 1e-04, ]$snp
#   ind.z <- z.gwas[[i]][z.gwas[[i]]$p.value < 1e-04, ]$snp
#   snps.use <- mclapply(unique(union(ind.y, ind.z)), ld.foo, LD) %>% unlist %>% unique
#
#   Gz <- Gz[, snps.use]
#   Gy <- Gy[, snps.use]
#
#   pcs <- pcs.prepare(Gz, Gy)
#   npcs <- min(which(pcs@vars > 0.8)[1], 6)
#   if(npcs < 2) npcs <- 2
#   thr <- pcs@vars[npcs - 1]
#
#   m1 <- pcs.model(pcs, group = 1, Y = w, family = "gaussian", threshold = thr)
#   m2 <- pcs.model(pcs, group = 2, Y = y, family = "gaussian", threshold = thr)
#   tmp <- coloc.test(m1, m2, bayes = TRUE)
#   chisqr <- tmp@result['chisquare']
#   pchisq(chisqr, tmp@result['n'] - 1, lower.tail = FALSE)
# }) %>% unlist
