##########
##Elastic Net on HapMap3
##########
##Data Treatment
##Here is your directory
(getwd())

##Change wd below to desired directory
setwd()

##Load up the expression table
##IF YOU ARE RUNNING ON WINDOWS YOU MAY NEED TO CHANGE '/' to '\' for getting directory
chb <- read.csv(file = '~/chb_p3.csv',header=TRUE)
##Fix the data
rownames(chb) <- chb$Normalization.REF
chb1 <- chb[,-1]
chb1 <-t(chb1)
##pull out the ids
ids <- data.frame( colnames(chb1) )
colnames(ids) <- c('x')

##Get the annotation information
genecodez <- read.csv('~/GeneSource.csv',na.strings = c('',NA))
##Filter the annotation information with only relevant data
z <- merge(ids,genecodez,by.x='x',by.y='X1Reporter')
zz <- z[,1:9]
sum(grepl('chrX',zz$Entry.entrez.))
##remove the x chromosome annotations
ind <- which(grepl('chrX',zz$Entry.entrez.))
zzz <- zz[-ind,]
##remove the empty annotations
z4 <- zzz[!is.na(zzz$Entry.entrez.),]
genecode <- z4
realids <- data.frame(z4$x)

##fix original expression data with the genes we have annotated data for
Rchb <- merge(realids,chb,by.x = 'z4.x',by.y = 'Normalization.REF')
##final exp data matrix
t.expdata <- t(Rchb)
expsamplelist <- rownames(t.expdata)


##read in the bim & fim & fam files
bim <- read.table('~/hapmap3_r2_b36_fwd.consensus.qc.poly.bim', header = FALSE)
newfam  <- read.table('~/hapmap3_r2_b36_fwd.consensus.qc.poly.fam',header = FALSE)
##Select for the column of ids
fam <- newfam$V2
samplelist <- intersect(fam,expsamplelist)
##get expression of samples with genotypes###
exp.w.geno <- t.expdata[samplelist,] 
explist <- colnames(exp.w.geno)


##Left out for cross hacking due to intensity and length of loading the file, caused my R to crash, 3.2 gb table to read, >1mil rows/cols ...Shyam
##out <- scan('~/out',what = 'character')
##X <- out[,samplelist]
## X <- t(X)

workingbest <- "working_10-foldCV_elasticNet_alpha.txt"
workingweight <- "_elasticNet_alpha_weights_.txt"

##########
##Statistical Analysis:
install.packages('glmnet',dependencies = TRUE)
library(glmnet)
k = 10 #number of folds for cross val
alpha = 0.5 #.5 alpha needed to run elastic net, 0 is lasso, and 1 is ridge regression

resultsarray <- array(0,c(length(explist),8))
dimnames(resultsarray)[[1]] <- explist
resultscol <- c("gene","alpha","cvm","lambda.iteration","lambda.min","n.snps","R2","pval")
dimnames(resultsarray)[[2]] <- resultscol
weightcol = c("gene","SNP","refAllele","effectAllele","beta")
write(resultscol,file=workingbest,ncolumns=8,sep="\t")
write(weightcol,file=workingweight,ncol=5,sep="\t")

set.seed(321)

##############################CANT RUN PAST THIS PAST THIS POINT WITHOUT HAVING UPLOADED OUT##############################################################

for(i in 1:length(explist)){
  cat(i,"/",length(explist),"\n")
  gene <- explist[i]
  geneinfo <- gencode[gene,]
  chr <- geneinfo[1]
  c <- substr(chr$V1,4,5)
  start <- geneinfo$V3 - 1e6 ### 1Mb lower bound for cis-eQTLS
  end <- geneinfo$V4 + 1e6 ### 1Mb upper bound for cis-eQTLs
  chrsnps <- subset(bim,bim[,1]==c) ### pull snps on same chr
  cissnps <- subset(chrsnps,chrsnps[,4]>=start & chrsnps[,4]<=end) ### pull cis-SNP info
  cisgenos <- X[,intersect(colnames(X),cissnps[,2])] ### pull cis-SNP genotypes
  if(is.null(dim(cisgenos))){
    bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
  }else{
    minorsnps <- subset(colMeans(cisgenos), colMeans(cisgenos,na.rm=TRUE)>0) ###pull snps with at least 1 minor allele###
    minorsnps <- names(minorsnps)
    cisgenos <- cisgenos[,minorsnps]
    if(is.null(dim(cisgenos)) | dim(cisgenos)[2] == 0){###effectively skips genes with <2 cis-SNPs
      bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
    }else{
      
      exppheno <- exp.w.geno[,gene] ### pull expression data for gene
      exppheno <- scale(exppheno, center=T, scale=T)  ###need to scale for fastLmPure to work properly
      exppheno[is.na(exppheno)] <- 0
      rownames(exppheno) <- rownames(exp.w.geno)
      
      ##run Cross-Validation over alphalist
      fit <- cv.glmnet(cisgenos,exppheno,nfolds=k,alpha=alpha,keep=T,parallel=F) ##parallel=T is slower on tarbell, not sure why
      
      fit.df <- data.frame(fit$cvm,fit$lambda,1:length(fit$cvm)) ##pull info to find best lambda
      best.lam <- fit.df[which.min(fit.df[,1]),] # needs to be min or max depending on cv measure (MSE min, AUC max, ...)
      cvm.best = best.lam[,1]
      lambda.best = best.lam[,2]
      nrow.best = best.lam[,3] ##position of best lambda in cv.glmnet output
      
      ret <- as.data.frame(fit$glmnet.fit$beta[,nrow.best]) # get betas from best lambda
      ret[ret == 0.0] <- NA
      bestbetas = as.vector(ret[which(!is.na(ret)),]) # vector of non-zero betas
      names(bestbetas) = rownames(ret)[which(!is.na(ret))]
      
      pred.mat <- fit$fit.preval[,nrow.best] # pull out predictions at best lambda
      
    }
  }
  if(length(bestbetas) > 0){
    res <- summary(lm(exppheno~pred.mat))
    genename <- as.character(gencode[gene,6])
    rsq <- res$r.squared
    pval <- res$coef[2,4]
    
    resultsarray[gene,] <- c(genename, alpha, cvm.best, nrow.best, lambda.best, length(bestbetas), rsq, pval)
    
    
    ### output best shrunken betas for PrediXcan
    bestbetalist <- names(bestbetas)
    bestbetainfo <- bim[bestbetalist,]
    betatable<-as.matrix(cbind(bestbetainfo,bestbetas))
    betafile<-cbind(genename,betatable[,2],betatable[,5],betatable[,6],betatable[,7]) ##output "gene","SNP","refAllele","effectAllele","beta"
    write(t(betafile),file=workingweight,ncolumns=5,append=T,sep="\t") # t() necessary for correct output from write() function
    
  }else{
    genename <- as.character(gencode[gene,6])
    resultsarray[gene,1] <- genename
    resultsarray[gene,2:8] <- c(NA,NA,NA,NA,0,NA,NA)
    
  }
  write(resultsarray[gene,],file=workingbest,ncolumns=8,append=T,sep="\t")
}


write.table(resultsarray,file="10-foldCV_elasticNet_alpha0.5_results.txt",quote=F,row.names=F,sep="\t")

