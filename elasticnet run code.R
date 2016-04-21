setwd('/home/comp383-group1/PredixcanData/')
expList = c('GIH_p3','JPT_p3','LWK_p3','MEX_p3','MKK_p3','YRI_p3')
chrNum = c(1:8,10:22)
outerBool = TRUE
innerBool = TRUE
save(list=c('expList','chrNum','outerBool','innerBool'),file="myobjects")
rm(list=ls())
load(file="myobjects")
while (outerBool) {
  
  
  while (innerBool) {
    expName = expList[1]
    
    k =10
    alpha = 0.5
    Q = chrNum[1]
    "%&%" = function(a,b) paste(a,b,sep="")
    library(glmnet)
    library(data.table)
    ##Load up the expression table
    ##IF YOU ARE RUNNING ON WINDOWS YOU MAY NEED TO CHANGE '/' to '\' for getting directory
    expNametxt = 'expression Data/' %&% expName %&% '_expression.txt'
    chrName = 'chr'%&% Q
    chb = data.frame(fread(input = expNametxt))
    ##Fix the data
    rownames(chb) <- chb$Normalization.REF
    chb1 <- chb[-1,-1]
    chb1 <-t(chb1)
    ##pull out the ids
    ids <- data.frame( colnames(chb1) )
    colnames(ids) <- c('x')
    t.expdata = chb1
    
    g = read.csv('output.csv')
    g$id = gsub('\\+AF8-','_',g$id)
    genecode = data.frame(g$chromosome,g$X.3,g$start,g$stop,g$ensembl,g$id)
    rownames(genecode) <- genecode[,6]
    genecode <- genecode[genecode[,1]==chrName,] ##pull genes on chr of interest
    t.expdata <- t.expdata[,intersect(colnames(t.expdata),rownames(genecode))] ###pull gene expression data w/gene info
    gencode = genecode
    colnames(gencode) = c('V1','V2','V3','V4','V5','V6')
    expsamplelist <- rownames(t.expdata) ###samples with exp data###
    
    bvm <- data.frame(fread('hapmap3_r2_b36_fwd.consensus.qc.poly.bim'))
    rownames(bvm) <- bvm$V2
    bim <- bvm[bvm[,1]==Q,]
    
    fvm  <- data.frame(fread('zip fils/hapmap3_r2_b36_fwd.consensus.qc.poly.fam',header = FALSE))
    
    fam = fvm$V2
    samplelist <- intersect(fam,expsamplelist)
    
    exp.w.geno <- t.expdata[samplelist,] ###get expression of samples with genotypes###
    explist <- colnames(exp.w.geno)
  
    gtxfile = '/home/comp383-group1/PredixcanData/Chromosomes/Chr' %&% Q %&% '.raw'
    GT = data.frame(fread(gtxfile))
    gtx = GT[,-(1:6)] 
    rownames(gtx) <- fam
    X <- gtx[samplelist,]
    na.change <- function(COLX){
      COLX[is.na(COLX)] <- mean(COLX,na.rm = TRUE)
      COLX
    }
    X <- apply(X,2,na.change)
    X<- data.frame(X)
    TEST<-gsub(pattern = '.*_r',replacement = 'r',x = colnames(X),perl = TRUE)
    TEST<-gsub(pattern = '_.*',replacement = '',x = TEST,perl = TRUE)
    (head(TEST))
    colnames(X) = TEST
    
    resultsarray <- array(0,c(length(explist),8))
    dimnames(resultsarray)[[1]] <- explist
    resultscol <- c("gene","alpha","cvm","lambda.iteration","lambda.min","n.snps","R2","pval")
    dimnames(resultsarray)[[2]] <- resultscol
    workingbest <- 'Results/working_' %&% expName %&% "_exp_" %&% 10 %&% "-foldCV_elasticNet_alpha" %&% 0.5 %&% "_" %&% 'hapmap' %&% chrName %&% ".txt"
    write(resultscol,file=workingbest,ncolumns=8,sep="\t")
    
    weightcol = c("gene","SNP","refAllele","effectAllele","beta")
    workingweight <- 'Results/working_' %&% expName  %&% "_elasticNet_alpha" %&% '0.5' %&% "_" %&% 'hapmap' %&% '_' %&% chrName %&% ".txt"
    write(weightcol,file=workingweight,ncol=5,sep="\t")
    
    set.seed(312)
    i = 27
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
        minorsnps <- subset(Matrix::colMeans(cisgenos), Matrix::colMeans(cisgenos,na.rm=TRUE) > 0) ###pull snps with at least 1 minor allele###
        minorsnps <- names(minorsnps)
        cisgenos <- cisgenos[,minorsnps]
        if(is.null(dim(cisgenos)) | dim(cisgenos)[2] == 0){###effectively skips genes with <2 cis-SNPs
          bestbetas <- data.frame() ###effectively skips genes with <2 cis-SNPs
        }else{
          
          exppheno <- exp.w.geno[,gene] ### pull expression data for gene
          exppheno <- scale(as.numeric(exppheno), center=T, scale=T)  ###need to scale for fastLmPure to work properly
          #exppheno <- scale(exppheno, center=T, scale=T)  ###need to scale for fastLmPure to work properly
          
          exppheno[is.na(exppheno)] <- 0
          
          rownames(exppheno) <- rownames(exp.w.geno)
          
          ##run Cross-Validation over alphalist
          fit <- cv.glmnet(as.matrix(cisgenos),as.matrix(exppheno),nfolds=k,alpha=alpha,keep=T,parallel=F) ##parallel=T is slower on tarbell, not sure why
          
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
    write.table(resultsarray,file= '/home/comp383-group1/PredixcanData/Results/' %&%expName %&% k %&% "-foldCV_elasticNet_alpha" %&% alpha %&% "_" %&% 'hapmap' %&% chrName %&% '.txt',quote=F,row.names=F,sep="\t")
      rm(list=ls())
    load(file="myobjects")
    chrNum = chrNum[-1]
    if (length(chrNum) == 0) {
      innerBool = FALSE
    }
    save(list=c('expList','chrNum','outerBool','innerBool'),file="myobjects")
  }
  
  rm(list=ls())
  load(file="myobjects")
  chrNum = c(1:22)
  expList = expList[-1]
  innerBool = TRUE
  if (length(expList) == 0) {
    outerBool = FALSE
    }
  save(list=c('expList','chrNum','outerBool','innerBool'),file="myobjects")
  load(file="myobjects")
}


