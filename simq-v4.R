source("~/Projects/coloc-cond-mask/load-hapdata.R")

## now we add 2 cvs

library(randomFunctions)
library(magrittr)
devtools::load_all("~/RP/coloc")
args <- getArgs(defaults=list(N=400,NSIM=2,NSNP=500,ld="highld",pop="eur",
e1="A", # A, AB, -
e2="same", # same, weaker, diff, -
t="A"), # A, AB, BA, AC, C
                numeric=c("N","NSIM","NCV"))
print(args)

## hdata <- readRDS(paste0("~/scratch/",args$ld,".RDS"))
GG <- readRDS("/home/cew54/share/Projects/twas/sample_geno.rds")

makeD <- function(y,X,best=NULL) {
    if(is.null(best)) {
        m <- snp.rhs.estimates(y~1,snp.data=X,family="Gaussian")
    } else {
        df <- data.frame(y=y,best=as(X[,best],"numeric"),row.names=rownames(X))
        m <- snp.rhs.estimates(y ~ .,data=df,snp.data=X,family="Gaussian")
    }
    nulls <- sapply(m,is.null)
    if(any(nulls)) 
        m <- m[!nulls]
    b <- sapply(m,"[[","beta")
    v <- sapply(m,"[[","Var.beta")
    z <- b/sqrt(v)
    list(N=args$N,
         MAF=dfsnps$maf[!nulls],
         beta=b,
         varbeta=v,
         type="quant",
         sdY=sd(y),
         snp=dfsnps$id[!nulls])
}
getmaxz <- function(D) {
    wh <- which.max(abs(D$beta)/sqrt(D$varbeta))
    D$snp[ wh ]
}
valmaxz <- function(D) {
    z <- abs(D$beta)/sqrt(D$varbeta)
    max(z)
}

amax <- function(x) {
    x[which.max(abs(x))]
}

simone <- function() {
## make genotypes
G <- GG[[sample(1:length(GG),1)]] # sample one genotype scenario
nr <- nrow(G)
    G1 <- round(G[sample(1:nr,args$N,replace=TRUE),]/2 + G[sample(1:nr,args$N,replace=TRUE),]/2)
    G2 <- round(G[sample(1:nr,args$N,replace=TRUE),]/2 + G[sample(1:nr,args$N,replace=TRUE),]/2)
    ## colnames(G1) <- colnames(G2) <- dfsnps$id
    rownames(G1) <- rownames(G2) <- paste0("I",1:args$N)
    X1 <- new("SnpMatrix",G1+1)
    X2 <- new("SnpMatrix",G2+1)

## causal variants
ABC=sample(1:ncol(G), 4); names(ABC) <- c("A","B","C", "-")
## effects at those variants, expression
beta <- replicate(3,amax(rnorm(10,0,0.15)))   %>% c(., 0) # beta_A, beta_B, beta_C
names(beta) <- names(ABC)
    
## outcomes - trait
cv <- ABC[ unlist(strsplit(args$t,"")) ]
y  <-  rnorm(args$N) + G1[,cv,drop=FALSE] %*% beta[names(cv)]
## outcomes - expr x 5
cv <- ABC[ unlist(strsplit(args$e1,"")) ]
z1  <-  rnorm(args$N) + G2[,cv,drop=FALSE] %*% beta[names(cv)]
py <- single.snp.tests(y,snp.data=X1)  %>% chi.squared(.,1)  %>% pchisq(.,df=1,lower=FALSE)
p1 <- single.snp.tests(z1,snp.data=X2)  %>% chi.squared(.,1)  %>% pchisq(.,df=1,lower=FALSE)
if(min(py,na.rm=TRUE)>1e-4 || min(p1,na.rm=TRUE)>1e-4)
return(NULL)

b4 <- beta
if(args$e2=="weaker")
b4 <- beta * runif(length(beta))
if(args$e2=="-")
cv <- ABC["-"]
if(args$e2=="diff")
cv <- ABC[ setdiff(names(ABC), unlist(strsplit(c(args$t,args$e1),""))) ]
z4 <- replicate(4,rnorm(args$N) + G2[,cv,drop=FALSE] %*% b4[names(cv)], simplify=FALSE)   %>% do.call("cbind",.)
z <- cbind(z1,z4)
LD <- cor(G[,ABC]); dimnames(LD) <- list(names(ABC),names(ABC))
psnps <- which(py<1e-4 | p1 < 1e-4)
maxr2 <- cor(G)[,psnps,drop=FALSE]^2  %>% apply(.,1,max)
pocket <- which(maxr2 > 0.2)
return(list(key=c(args, list(e4=paste(names(cv),collapse=""),
beta=beta,b4=b4,LD=LD,ABC=ABC)),
                data=list(y=y,Gy=X1,z=z,Gy=X2,pocket=pocket)))
} 

for(i in 1:args$NSIM) {
x <- NULL
while(is.null(x))
    x=simone()
    ## dir.create("~/share/Projects/twas/sims",recursive=TRUE)
    ## dir.create("~/scratch/Projects/twas/sims",recursive=TRUE)

    tmp <- tempfile(pattern="simv4","~/share/Projects/twas/sims",fileext=".rds")
    key <- sub("sims","keys",tmp)
    
    saveRDS(x$data, file=tmp, version=2)
    saveRDS(x$key, file=key, version=2)
}


