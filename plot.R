## source("~/Projects/coloc-cond-mask/load-hapdata.R")
setwd("~/D")
library(randomFunctions)
library(magrittr)
library(data.table)
library(progress)
library(viridis)
library(ggplot2)
library(cowplot); theme_set(theme_cowplot())
## devtools::load_all("~/RP/coloc")
library(network)
library(sna)
library(ggnet)

res <- readRDS(file = '~/share/Projects/twas/sim_results-v5.rds')
head(res)
## res[,ind.effect:=ifelse(effect!=0, "eqtl effect", "no effect")]
## res[,match:=ifelse(cv=="same","effects same","effects diff")]
## res[,summary(coloc.pval),by="version"]
with(res,table(t,e1))
with(res,table(t,e4))

res[,pname:=paste(t,e1,sub("-","",e4),sep=":")]
## E1 is the test trait
## trait = E1-E5 are the 5 expression traits tested
## t is the causal variant for the potentially colocalising trait
## e1 is the causal variant for E1
## e4 is the causal variant for the other E
## e2 is an indicator - are they the same between E1 and E2-E5, same but weaker, or missing (no eqtl) or different
## match is whether t and E1 match causal variants
## background is whether E1 and E2-E5 match causal variants

## networks

## draw this:
library(hrbrthemes)
cols <- ipsum_pal()(9)  %>% c(.,"#aaaaaa","#ffffff")
library(igraph)
patterns <- unique(res[,.(t,e1=sub("-","",e1),e4=sub("-","",e4),pname)])
nodes <- data.table(name=c("GWAS","A","B","C","D","E.test","E.back"),
                    class=c("GWAS","v","v","v","v","expression","expression"),
                    col=cols[c(2,10,10,10,10,2,4)],
                    x=c(1.5,1,2,3,4,3.5,2.5),
                    y=c(3,2,2,2,2,3,3),
                    stringsAsFactors=FALSE)
layout <- matrix(c(nodes$x,nodes$y),ncol=2,dimnames=list(nodes$name,c("x","y")))
makerel <- function(p,to) {
  rels <- lapply(p, function(x)
    data.frame(from=if(x=="") {
                      ""
                    } else {
                      unlist(strsplit(x,""))
                    },
               to=to))
  names(rels) <- p
  rels
}
  trel <- makerel(unique(patterns$t),to="GWAS")
e1rel <- makerel(unique(patterns$e1),to="E.test")
e4rel <- makerel(unique(patterns$e4)  %>%  setdiff(.,""),to="E.back")

## rows <- data.frame(from=c("A","B"),to=c("B","C"))
## grows <- graph_from_data_frame(rows, directed=FALSE, vertices=nodes) 

##   ## use network, doesn't tile well
## plotter <- function(i) {
##   relations <- rbind(trel[[ patterns$t[i] ]],
##                      e1rel[[ patterns$e1[[i]] ]],
##                      e4rel[[ patterns$e4[[i]] ]])  %>% as.matrix()
##   adj <- matrix(0,7,7,dimnames=list(nodes$name,nodes$name))
##   adj[relations] <- 1
##   net <- network(adj)
##   net %v% "class" = nodes$class
##   net %v% "name" = nodes$name
##   pal <- nodes$col
##   names(pal) <- nodes$name
##   net %v% "color" = nodes$col
##   net %v% "x" = nodes$x
##   net %v% "y" = nodes$y
##   edges <- as.matrix(net,matrix.type="edgelist")
##   net %e% "color" = nodes$col[ edges[,2] ] # colour by node the edge is to
##   ggnet2(net,mode=c("x","y"),
##          color="color",edge.color="color",
##          ## arrow.size=10, arrow.gap=0.08,
##          label="name",
##          size=20,
##          shape=15) +
##     theme(axis.line=element_blank())
## }
## plotter(1)
## glist <- lapply(1:nrow(patterns), plotter)

  ## cook my own
plotter <- function(i) {
  ti <- "No match"
  if(patterns$t[i]==patterns$e1[i]) {
    ti <- "Full match"
  } else {
    if(length(intersect(strsplit(patterns$t[i],"")[[1]],
                        strsplit(patterns$e1[i],"")[[1]])))
      ti <- "Partial match"
  }
  relations <- rbind(trel[[ patterns$t[i] ]],
                     e1rel[[ patterns$e1[[i]] ]],
                     e4rel[[ patterns$e4[[i]] ]])  %>% as.matrix()
  rdf <- merge(relations,nodes[,.(name,x,y)],by.x="from",by.y="name")
  rdf <- merge(rdf,nodes[,.(name,x,y,col)],by.x="to",by.y="name",suffixes=c(".from",".to"))
  cscale <- structure(unique(nodes$col),names=unique(nodes$col))
  ggplot(nodes, aes(x=x,y=y)) +
    geom_segment(aes(x=x.from,y=y.from,xend=x.to,yend=y.to,col=col),data=rdf,size=1) +
    geom_point(aes(colour=col,size=class),pch=20) +
    geom_text(aes(label=name)) +
    xlim(0,5) + ylim(1.9,3.5) +
    scale_colour_manual(values=cscale) +
    scale_size_manual(values=c(GWAS=25,expression=25,v=10)) +
    ggtitle(ti) +
    theme(legend.position="none",axis.line=element_blank(),
          axis.title=element_blank(),
          axis.ticks=element_blank(), axis.text=element_blank())
} 
plotter(1)

  glist <- lapply(1:nrow(patterns), plotter)
    
## plot results by pattern
res[,lasso.twas.all:=ifelse(is.na(lasso.twas),1,lasso.twas)]
setnames(res,"lasso.twas","lasso.twas.tested")
mvars <- grep("twas",names(res),value=TRUE)
m <- melt(res[trait %in% c("E1")],
          id.vars = c("trait","t","e1","e4","coloc.pval","pname"),
          measure.vars=mvars)
## m <- m[!is.na(value)]
## m[!is.na(value),fdr:=p.adjust(value,method="BH"),by=c("variable","trait")]
## m[is.na(value),value:=1]
## m[is.na(fdr),fdr:=1]
head(m)

myCols <- viridis(option = "D", 4)

## equivalent:
## AB:A == AC:A
m[pname=="AC:A",pname:="AB:A"]

ggplot(res, aes(x=-log10(min.pval.expr), fill=!is.na(lasso.twas.tested))) +
  geom_histogram(position="dodge") +
  scale_fill_discrete("lasso run") +
  facet_wrap(~e1)

  geom_violin() + facet_wrap(~pname) + theme(legend.position="none")

ggplot(res, aes(x=is.na(lasso.twas), y=-log10(min.pval.expr), fill=is.na(lasso.twas))) +
  geom_violin() + facet_wrap(~pname) + theme(legend.position="none")

## summarised
ms <- m[,list(avg=mean(value<0.05,na.rm=TRUE)),by=c("t","e1","e4","pname","variable")]
ms[,variable:=sub(".twas","",variable)]
nrow(ms)

plotbar <- function(p) {
  msub <- ms[pname==p]
  ggplot(msub) +
    geom_col(aes(x=variable,y=avg,fill=variable)) +
    geom_hline(yintercept=0.05) +
    ylim(0,max(ms$avg,na.rm=TRUE)) +
    background_grid(major="y") +
    theme(legend.position="none",axis.line=element_blank())
}
bars <- lapply(patterns$pname, plotbar)

allplots <- lapply(1:nrow(patterns), function(i)
  plot_grid(glist[[i]], bars[[i]],ncol=1,rel_heights=c(1,2)))
plot_grid(plotlist=allplots)
  
plotg <- function(g) {
  x <- as.data.frame(g)
  
library(cowplot)
p3 <- ~plot(mtcars$qsec, mtcars$mpg)

plot_grid(~glist[[1]], bars[[1]], labels = c('A', 'B'), label_size = 12)

p
tab <- msub$avg
names(tab

ggplot(ms) +
  geom_col(aes(x=variable,y=avg,fill=variable)) +
  geom_hline(yintercept=0.05) +
  facet_wrap(~pname)

ggplot(ms[],aes(x=fdr_bin,y=avg,fill=variable)) +
    ## geom_histogram(position="dodge",binwidth=0.1,center=0.05) +
  geom_col(position="dodge") +
  ## facet_grid(match~background,labeller=label_both, scales="free") +
  facet_wrap(~pname,scales="free") +
  scale_fill_manual(values = myCols) +
    ggtitle("twas FDR")


                                        #DENSITY plot
ggplot(m[],aes(x=fdr,fill=variable)) +
  geom_histogram(aes(y = ..density..), position="dodge",binwidth=0.1,center=0.05) +
  ## facet_grid(match~background,labeller=label_both, scales="free") +
  facet_wrap(~pname,scales="free") +
  scale_fill_manual(values = myCols) +
  ggtitle("twas FDR")

#STANDARD plot
ggplot(m[],aes(x=fdr,fill=variable)) +
    geom_histogram(position="dodge",binwidth=0.1,center=0.05) +
  ## facet_grid(match~background,labeller=label_both, scales="free") +
  facet_wrap(~pname,scales="free") +
  scale_fill_manual(values = myCols) +
    ggtitle("twas FDR")


    

#MATCH
with(res[trait=="E1"],table(e1,t))
res[,match:="diff"]
res[trait=="E1" & e1=="-", match:="none"]
res[trait!="E1" & e4=="-", match:="none"]
res[trait=="E1" & t!="C" & e1!="-", match:="partial"]
res[trait!="E1" & e4 == "A" & t %in% c("AB", "AC"), match := "partial"]
res[trait!="E1" & e4 == "AB" & t %in% c("A", "AC"), match := "partial"]
res[trait!="E1" & e4 == "B-" & t == "AB", match := "partial"]
res[trait!="E1" & e4 == "BC-" & t %in% c("AB", "AC", "C"), match := "partial"]
res[trait=="E1" & t==e1, match:="full"]
res[trait!="E1" & t==e4, match:="full"]

with(res[trait == "E1"], table(e1, t, match))
with(res[trait != "E1"], table(e4, t, match))

#BACKGROUND
res[,background:="diff"]
res[e1 == '-' | e4 == '-', background:="none"]
res[e4 == e1 & background != "none", background:="full"]
res[e1 == "A" & e4 == "AB", background := "partial"]
res[e1 == "AB" & e4 %in% c("B-", "BC-"), background := "partial"]
with(res, table(e1, e4, background))

with(res[trait %in% c("E1", "E2")], table(background, match))


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

#--------------------->>>>>>PLOTS<<<<<<--------------------------------
myCols <- viridis(option = "D", 4)
#DENSITY plot
ggplot(m[keep==TRUE],aes(x=fdr,fill=variable)) +
  geom_histogram(aes(y = ..density..), position="dodge",binwidth=0.1,center=0.05) +
  facet_grid(match~background,labeller=label_both, scales="free") + scale_fill_manual(values = myCols) +
  ggtitle("twas FDR")

#STANDARD plot
ggplot(m[keep==TRUE],aes(x=fdr,fill=variable)) +
    geom_histogram(position="dodge",binwidth=0.1,center=0.05) +
    facet_grid(match~background,labeller=label_both, scales="free") + scale_fill_manual(values = myCols) +
    ggtitle("twas FDR")

    
