## source("~/Projects/coloc-cond-mask/load-hapdata.R")

library(randomFunctions)
library(magrittr)
library(data.table)
## devtools::load_all("~/RP/coloc")

dkeys <- "~/share/Projects/twas/keys"
ddata <- sub("keys","sims",dkeys)  %>% paste0(.,"/results")
keyfiles <- list.files(dkeys,pattern="^simv4")
datafiles <- list.files(ddata,pattern="^simv4")

files <- intersect(keyfiles,datafiles)
message("keyfiles found: ",length(keyfiles))
message("matching output files found: ",length(files))

files <- sample(files,200)
res <- vector("list",length(files))

library(progress)

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
head(res)
## res[,ind.effect:=ifelse(effect!=0, "eqtl effect", "no effect")]
## res[,match:=ifelse(cv=="same","effects same","effects diff")]
## res[,summary(coloc.pval),by="version"]

with(res,table(t,e1))
with(res,table(t,e4))

## E1 is the test trait
## trait = E1-E5 are the 5 expression traits tested
## t is the causal variant for the potentially colocalising trait
## e1 is the causal variant for E1
## e4 is the causal variant for the other E
## e2 is an indicator - are they the same between E1 and E2-E5, same but weaker, or missing (no eqtl) or different
## match is whether t and E1 match causal variants
## background is whether E1 and E2-E5 match causal variants

with(res[trait=="E1"],table(e1,t))
res[,match:="diff"]
res[trait=="E1" & e1=="-", match:="none"]
res[trait!="E1" & e4=="-", match:="none"]
res[trait=="E1" & t!="C" & e1!="-", match:="partial"]
res[trait!="E1" & t!="C" & e4!="-", match:="partial"]
res[trait=="E1" & t==e1, match:="full"]
res[trait!="E1" & t==e4, match:="full"]

with(res[trait=="E1"],table(e1,e4))
res[,background:="diff"]
res[trait!="E1" & e1=="-", background:="none"]
res[trait=="E1" & e4=="-", background:="none"]
res[trait!="E1" & e1!="-", background:="partial"]
res[trait=="E1" & e4!="-", background:="partial"]
res[trait!="E1" & e4==e1, background:="full"]
res[trait=="E1" & e1==e4, background:="full"]

with(res[trait %in% c("E1","E2")],table(background,match))

with(res,table(e1,e4,t))


mvars <- grep("twas",names(res),value=TRUE)
m <- melt(res[trait %in% c("E1","E2")],
          id.vars = c("trait","match","background","coloc.pval"),
          measure.vars=mvars)
head(m)
## m[,val0:=ifelse(is.na(value),runif(nrow(m)),value)]
m[,fdr:=p.adjust(value,method="BH"),by=c("variable","trait")]
## m <- m[version=="v3"]

tt <- with(m, table(match, background)) # max is 75000, min is 2358
tt

## m[,n:=.N,by=c("match","background","variable")]
## m[,keep:=sample(1:.N)<=max(tt),by=c("match","background")]
m[,keep:=TRUE]


library(ggplot2)
## m[,cond:=sub("all","equal strength",cond)]
## ggplot(m,aes(x=val0,fill=match)) + geom_histogram(position="dodge",binwidth=0.1) + facet_grid(variable~ind.effect+cond) +
##   ggtitle("twas p value, NA set to runif()")

ggplot(m[keep==TRUE],aes(x=value,fill=variable)) +
  geom_histogram(position="dodge",binwidth=0.1,center=0.05) +
  facet_grid(match~background,labeller=label_both, scales="free") +
  ggtitle("twas p values")

ggplot(m[keep==TRUE],aes(x=fdr,fill=variable)) +
  geom_histogram(position="dodge",binwidth=0.1,center=0.05) +
  facet_grid(match~background,labeller=label_both, scales="free") +
  ggtitle("twas FDR")


## ggplot(m[keep==TRUE],aes(x=fdr,fill=match)) +
##   geom_histogram(position="dodge",binwidth=0.1) +
##   facet_grid(variable~ind.effect+cond) +
##   ggtitle("twas FDR")

## ggplot(m0, aes(x=coloc.pval,fill=match)) +
##   geom_histogram() +
##   facet_grid(variable ~ ind.effect + cond) +
##   ggtitle("coloc p value distribution according to effects same or diff",
## subtitle="amongst FDR<0.01 TWAS")


## ggplot(m0, aes(x=coloc.pval,fill=match)) +
##   geom_histogram(aes(x=coloc.pval,stat(density))) +
##   facet_grid(variable ~ ind.effect + cond) +
##   ggtitle("coloc p value distribution according to effects same or diff",
## subtitle="amongst FDR<0.01 TWAS")



## ggplot(m,aes(x=variable,y=-log10(val0))) + geom_boxplot() + facet_grid(cond~ind.effect)

res[,mtl.fdr:=p.adjust(mtl.twas,method="BH")]
res[,fuser.fdr:=p.adjust(fuser.twas,method="BH")]
res[,rf.fdr:=p.adjust(rf.twas,method="BH")]
res[,lasso.fdr:=p.adjust(lasso.twas,method="BH")]

res[,fdr.sig:=pmin(mtl.fdr,fuser.fdr,rf.fdr,lasso.fdr,na.rm=TRUE)<0.01]
## ggplot(res[ind.effect=="eqtl effect"], aes(x=coloc.pval)) +
##   geom_histogram(binwidth=0.05) +
##   facet_grid(fdr.sig~match) +
##   ggtitle("coloc p values")


## ggplot(m, aes(x=coloc.pval,fill=match)) + geom_histogram() +
##   facet_grid(variable ~ ind.effect + cond) +
##   ggtitle("coloc p value distribution according to effects same or diff",
## subtitle="amongst all tests")


m0 <- m[fdr<0.01][order(coloc.pval)]
m0[,N:=1:.N,by=c("variable","trait")]
m0[,Ndiff:=cumsum(match=="diff"),by=c("variable","trait","ind.effect","cond")]
m0[,coloc.fdr:=(N-Ndiff)/N]

ggplot(m0, aes(x=coloc.pval,fill=match)) + geom_histogram() +
  facet_grid(variable ~ ind.effect + cond) +
  ggtitle("coloc p value distribution according to effects same or diff",
subtitle="amongst FDR<0.01 TWAS")



ggplot(m0,aes(x=coloc.pval,y=coloc.fdr)) +
  geom_path() +
  geom_smooth(se=FALSE) +
  facet_grid(variable ~ ind.effect+cond) +
  ggtitle("coloc false acceptance of same causal variants")

res[,N:=cumsum(1:.N),

alpha <- seq(0,0.9,by=0.01)
N <- sapply(alpha, function(a) 
    nrow(res[fdr.sig==TRUE & coloc.pval>a]))
P <- res[fdr.sig==TRUE]$coloc.pval
FDR <- sapply(alpha, function(a) 
    with(res[fdr.sig==TRUE & coloc.pval>a], mean(match=="effects diff")))
FRR <- sapply(alpha, function(a) 
    with(res[fdr.sig==TRUE & coloc.pval>a], mean(match=="effects same")))
par(mfrow=c(1,2))
hist(P)
#plot(alpha,N * (1-FDR))
plot(alpha,FDR)

summary(N)

plot(alpha,FRR)

with(res[fdr.sig==TRUE ], mean(match=="effects diff"))
