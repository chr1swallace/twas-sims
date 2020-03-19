source("~/Projects/coloc-cond-mask/load-hapdata.R")

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

library(mvtnorm)

simone <- function() {
  ## make genotypes
  G <- GG[[sample(1:length(GG),1)]] # sample one genotype scenario
  nr <- nrow(G)
  ## these have similar correlation structures and slightly more
  ## variable MAF than the input data, but constructed in same way for
  ## expression and GWAS, so valid hack.
  G1 <- round(G[sample(1:nr,args$N,replace=TRUE),]/2 +
              G[sample(1:nr,args$N,replace=TRUE),]/2)
  G2 <- round(G[sample(1:nr,args$N,replace=TRUE),]/2 +
              G[sample(1:nr,args$N,replace=TRUE),]/2)
  drop <- apply(G1,2,var)==0 | apply(G2,2,var)==0 # just in case
  G1 <- G1[,!drop]
  G2 <- G2[,!drop]
  G <- G[,!drop]
  ## colnames(G1) <- colnames(G2) <- dfsnps$id
  rownames(G1) <- paste0("I",1:args$N)
  rownames(G2) <- paste0("I",1:args$N)
  X1 <- new("SnpMatrix",G1+1)
  X2 <- new("SnpMatrix",G2+1)

  ## causal variants
  ABC=sample(1:ncol(G), 5); names(ABC) <- c("A","B","C","D","-")
  ## effects at those variants, use amax() of 10 draws to sample from
  ## a bimodal distribution with lower pdf at 0, as these draws will
  ## be discarded by the "will we run this" filters below.
  ##
  ## Idea of sampled distribution
  ## > replicate(10000,amax(rnorm(10,0,0.15))) %>% abs() %>% summary()
  ##  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  ## 0.05418 0.22635 0.27188 0.28027 0.32649 0.65102
  beta <- replicate(4,amax(rnorm(5,0,0.2)))   %>% c(., 0) # beta_A, beta_B, beta_C
  names(beta) <- names(ABC)
    
  ## outcomes - expr 
  cv <- ABC[ unlist(strsplit(args$e1,"")) ]
  S <- matrix(args$rho,5,5)
  diag(S) <- 1
  Z <- rmvnorm(args$N,sigma=S)
  z1  <-  Z[,1] + G2[,cv,drop=FALSE] %*% beta[names(cv)]
  p1 <- single.snp.tests(z1,snp.data=X2)  %>% chi.squared(.,1)  %>% pchisq(.,df=1,lower=FALSE)

  ## outcomes - expr background
  ## add some wobble to beta
  b4 <- beta * runif(length(beta),min=0.8,max=1.2)
  if(args$e2=="-")
    cv <- ABC["-"]
  if(args$e2=="diff") {
    used <- c(args$t,args$e1)  %>% strsplit(.,"")  %>% unlist()  %>% unique()
    cv <- ABC[ setdiff(names(ABC), c(used,"-")) ][1]
  }
  if(args$e2=="t")
    cv <- ABC[ unlist(strsplit(args$t,"")) ]
  z4 <- Z[,-1] + matrix(G2[,cv,drop=FALSE] %*% b4[names(cv)],nrow=args$N,ncol=ncol(Z)-1)
e4 <- cv # store
  p2 <- single.snp.tests(z4[,1],snp.data=X2)  %>% chi.squared(.,1)  %>% pchisq(.,df=1,lower=FALSE)
  p3 <- single.snp.tests(z4[,2],snp.data=X2)  %>% chi.squared(.,1)  %>% pchisq(.,df=1,lower=FALSE)
  p4 <- single.snp.tests(z4[,3],snp.data=X2)  %>% chi.squared(.,1)  %>% pchisq(.,df=1,lower=FALSE)
  p5 <- single.snp.tests(z4[,4],snp.data=X2)  %>% chi.squared(.,1)  %>% pchisq(.,df=1,lower=FALSE)

  ## rule: only run if at least one expression p1..p4 < 1e-7
  minp <- min(c(p1,p2,p3,p4,p5),na.rm=TRUE)
  if(minp > 1e-7)
    return(NULL)
  z <- cbind(z1,z4)

  ## outcomes - trait
  cv <- ABC[ unlist(strsplit(args$t,"")) ]
  y  <-  rnorm(args$N) + G1[,cv,drop=FALSE] %*% beta[names(cv)]
  py <- single.snp.tests(y,snp.data=X1)  %>% chi.squared(.,1)  %>% pchisq(.,df=1,lower=FALSE)

  minp <- c(y=min(py),e1=min(p1),e2=min(c(p2,p3,p4,p5)))

  
  LD <- cor(G[,ABC]); dimnames(LD) <- list(names(ABC),names(ABC))
  psnps <- c(which(py<1e-4 | p1 < 1e-4))  %>% unique()
  if(length(psnps)) {
    maxr2 <- cor(G,G[,psnps,drop=FALSE])^2 %>% apply(.,1,max)
    pocket <- which(maxr2 > 0.2)
  } else {
    pocket <- numeric(0)
  }
  return(list(key=c(args, list(e4=paste(names(e4),collapse=""),
                               beta=beta,b4=b4,LD=LD,ABC=ABC,minp=minp)),
                data=list(y=y,Gy=X1,z=z,Gy=X2,pocket=pocket)))
} 

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


