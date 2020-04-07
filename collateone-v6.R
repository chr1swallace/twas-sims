library(data.table)
library(magrittr)
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())

setwd("~/D")
source("seaborn_pal.R")
files <- list.files(pattern="simv6",path="~/share/Projects/twas/simone",full=TRUE)
message("files found: ",length(files))
if(length(files)>1000 & interactive()) # for testing
  files <- sample(files,1000)
f=files[1]

if(!file.exists("~/collateone.RData")) {
  library(progress)
  pb <- progress_bar$new(total = length(files))
  results <- vector("list",length(files))
  for(i in seq_along(files)) {
    pb$tick()
    f <- files[i]
    data=readRDS(f)
    ## str(data)
    betas <- as.list(data$key$beta[ unlist(data$key[c("e1","t")]) ])
    names(betas) <- c("beta.e1","beta.t")
    results[[i]] <- c(data$key[c("e1","t")],
                      betas,
                      as.list(data$twas))
  }
  results %<>% rbindlist(.,fill=TRUE)
  save(results, file="~/collateone.RData")
} else {
  load(file="~/collateone.RData")
}

## lasso constant fit
head(results)
summary(results)
table(results$e1,!is.na(results$lasso.p))
results[,eQTL:=sub("-","none",e1)]
results[,eQTL:=sub("A","effect",eQTL)]
results[,eQTL:=sub("B","different",eQTL)]
p.constants <- ggplot(results[e1 %in% c("-","A")],aes(x=eQTL,fill=is.na(lasso.p))) +
  geom_bar(position="fill",col="grey") +
  scale_fill_seaborn("Constant prediction",palette="pastel",direction=-1) +
  background_grid(major="y") +
  labs(y="Proportion",x="eQTL causal variant") +
  theme(legend.position="right",legend.direction="vertical",axis.line.y=element_blank())
p.constants

## effect distribution
tmp <- results[e1=="A"]
## ## tmp[,truebeta:=beta.t/beta.e1]
## ggplot(tmp[lasso.p<1e-4],aes(x=lasso.est)) + geom_histogram() + geom_vline(xintercept=1,col="red")
## ggplot(tmp[rf.p<1],aes(x=rf.est)) + geom_histogram() + geom_vline(xintercept=1,col="red")
## ggplot(tmp,aes(x=rf.est,y=-log10(rf.p))) + geom_point() + geom_vline(xintercept=1,col="red")
m <- melt(tmp,measure.vars=list(est=c("lasso.est","rf.est"),
                                se=c("lasso.se","rf.se"),
                                p=c("lasso.p","rf.p")
                                ))
head(m)
m[,method:=c("Lasso","RF")[variable]]

## p.effects <- ggplot(m,aes(y=est,ymin=est-1.96*se,ymax=est+1.96*se,
##                x=-log10(p),col=method)) +
##   geom_point(alpha=0.2) +
##   geom_linerange(lwd=0.2,alpha=0.2) +
##   geom_hline(yintercept=1,col="black") + 
##   geom_smooth(se=FALSE) +
##   ## facet_grid(method~.) +
##   labs(y="Ratio of estimated to true TWAS effect") +
##   xlim(0,40) +
##   ## scale_colour_seaborn("Method") +
##   scale_colour_manual(values=myCol) +
##   background_grid(major="y") +
##   ## scale_y_log10() +
##   theme(legend.position="bottom",axis.line.y=element_blank())

use <- c(sample(which(m$est>1.96*m$se & m$method=="Lasso"),2000),
         sample(which(m$est>1.96*m$se & m$method=="RF"),2000))
p.effects.log <- ggplot(m[use],
                        aes(y=est,ymin=est-1.96*se,ymax=est+1.96*se,
               x=-log10(p),col=method)) +
  geom_point(alpha=0.2,size=0.2) +
  geom_linerange(lwd=0.1,alpha=0.2) +
  geom_hline(yintercept=1,col="black") + 
  geom_smooth(se=FALSE) +
  ## facet_grid(method~.) +
  labs(y="Relative TWAS effect estimate") +
  labs(x="-log10 P value",limits=c(0,40)) +
  ## scale_colour_seaborn("Method") +
  scale_colour_manual(values=myCol) +
  background_grid(major="y") +
  scale_y_log10() +
  theme(legend.position="bottom",axis.line.y=element_blank())


## p value distribution
m <- melt(results,measure.vars=c("rf.p","lasso.p"))
m[variable=="rf.p",variable:="RF"]
m[variable=="lasso.p",variable:="Lasso"]
m[,variable:=factor(as.character(variable))]
m[is.na(value),value:=1]

ggplot(m,aes(x=value,fill=variable)) +
  geom_histogram(position="dodge") +
  facet_grid(e1 ~ .,scales="free_y") +
  ## scale_fill_seaborn("Method") +
  scale_fill_manual(values=myCol) +
  labs(x="P value") +
  background_grid(major="y") +
  theme(legend.position="bottom")

p.pvals <- ggplot(m[e1=="A"],aes(x=-log10(value))) +
  geom_histogram(aes(fill=variable),position="dodge") +
  ## geom_histogram(aes(fill=variable),position="dodge",stat="density") +
  ## geom_density(aes(col=variable)) +
  ## facet_grid(e1 ~ .,scales="free_y") +
  labs(x="-log10 P value",limits=c(0,40)) +
  ## scale_fill_seaborn("Method") +
  ## scale_colour_seaborn("Method") +
  scale_colour_manual(values=myCol) +
  scale_fill_manual(values=myCol) +
  background_grid(major="y") +
  theme(legend.position="bottom",axis.line.y=element_blank())

top <- plot_grid(p.constants,p.pvals,labels=c("a","b"),nrow=1,rel_widths=c(.4,.6))
plot_grid(p.pvals,
          p.effects.log,
          nrow=2,#rel_widths=c(.2,.6,.2),
          labels=c("a","b"),
          align="x")
## Effects of lasso shrinkage on TWAS. a Lasso produces constant predictions when there is no effect ("none") as well as in a proportion of simulations where an effect is present. b Lasso-TWAS p values amongst simulations with shared eQTL/GWAS causal variants show a spike at p=1, and a longer tail than RF, indicating that weaker effects are missed, but that stronger effects show greater significance. c TWAS effect estimates (effect of expression on GWAS trait) are underestimated for weak effects for RF, tending to 1 for stronger effects. For lasso, TWAS effect estimates are systematically over estimated, even for well-powered studies.
ggsave("sims-lasso-rf.png",height=6,width=8,scale=1.3)

################################################################################

## junk below here
m[,.(p05=mean(value<0.05),
     p01=mean(value<0.01),
     p005=mean(value<0.005),
     p001=mean(value<0.001)),
    ,by=c("variable","e1")]

## fdr distribution
nulls <- results[e1=="-"]
ab <- results[e1!="-"]
a <- results[e1=="A"]
b <- results[e1=="B"]
n <- nrow(ab)

f <- 0.2
fakemix <- function(f) { # 
  fa <- rbind(nulls[ sample(1:nrow(nulls), round(n*(1-f)), replace=TRUE) ],
              ab[ sample(1:nrow(ab), round(n*(f))) ])
  m <- melt(fa,measure.vars=c("rf.p","lasso.p"))
  m[is.na(value),value:=1]
  m[,fdr:=p.adjust(value,method="BH"),by="variable"]
  m[e1=="A",.(f=f,
              fdr05=sum(fdr<0.05),
              fdr01=sum(fdr<0.05)),by="variable"]
}
fakes <- lapply(seq(0.001,0.1,by=0.001),fakemix)  %>% rbindlist()

ggplot(fakes,aes(x=f,y=fdr05,col=variable)) +
  geom_point() +
  geom_path() +
  scale_colour_seaborn() +
  background_grid(major="y") +
  theme(legend.position="bottom")

ggplot(fakes,aes(x=f,y=fdr01,col=variable)) +
  geom_point() +
  geom_path() +
  scale_colour_seaborn() +
  background_grid(major="y") +
  theme(legend.position="bottom")

m[,fdr:=p.adjust(value,method="BH")]
m[,mean(fdr<0.05),by=c("variable","e1")]
