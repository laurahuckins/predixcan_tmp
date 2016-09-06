library(glmnet)

library(parallel)
stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(x))
lower <- function(x) quantile(x,0.025,na.rm=TRUE)
upper <- function(x) quantile(x,0.975,na.rm=TRUE)


read.files<-function(IND, chr, code, dir){

    if(IND==1){ #  RNAseq file

        output<-read.delim("/sc/orga/projects/CommonMind/results/phaseI/DIFF_GENE_EXPRESSION/OFFICIAL_RUNS/ensembl.REQUIRE_ANCESTRY.NO_HIDDEN_CONFOUND.STANDARD_COVARS/DLPFC.ensembl.KNOWN.ADJUSTED.VOOM_NORMALIZED.GE.WEIGHTED_RESIDUALS.tsv", header=T, sep='')
        rownames(output)<-output[,1]
        header=read.delim(sprintf("%s/header.rnaseqfile",dir))
        colnames(output)<-toupper(colnames(header))

     } else if (IND == 2){ # snplist

        output=read.delim("/sc/orga/work/huckil01/filt_file/snplist_noindels.gz", header=F)
            
    } else if (IND == 3){ # List of genes for this chr only

        gencode.perm<-read.delim("/sc/orga/projects/CommonMind/lhuckins/prediXcan_files/gencode.v18.genes.patched_contigs.summary.protein.named.txt", header=F, sep='')
        genes.use<-gencode.perm[which(gencode.perm[,1]==chr),9]
        output<-genes.use

    } else if (IND == 4){ # Dosage data

        output<-read.delim(sprintf("/sc/orga/projects/CommonMind/lhuckins/CM5-chr%s_imputed.dos.id.gz",chr), header=T, sep='')

        dupl<-which(duplicated(output[,2])==TRUE)
        if(length(dupl)>0){
                 output<-output[-(dupl),]
        }
        rownames(output)<-output[,2]
        colnames(output)<-toupper(colnames(output))

    } else if (IND == 5){ # Cross-validation folds

        output=vector(10, mode="list")    # This is so I can calc for each set of indivs at the same time

        for (i in 1:10){
    
            output[[i]]<-read.delim(sprintf("/sc/orga/projects/CommonMind/lhuckins/prediXcan_files/BSLMM/Ids.set%s.use", i), header=F, sep='')

        }

    } else if (IND == 6){ # All snp lists for all genes

        output<-system("ls /hpc/users/huckil01/snp_lists/E*.onemb.gz | sed 's/.*E/E/g;s/.onemb.gz//g' ", intern=T)


    } else if (IND == 7){ # read bimfile, get only lines for this chr

        bimfile<-read.delim("/sc/orga/work/huckil01/filt_file/CM-pos-imputed-hardcall95.maf5.info08.ld08.nodupls.noindels.noambig.bim.gz", header=F, sep='')
        bimfile<-bimfile[which(bimfile[,1]==chr),]
        rownames(bimfile)<-bimfile[,2]
        output<-bimfile
    }

    return(output)
                                    
}


output<-mclapply(1:8, read.files, chr, code, dir)


rnaseqfile<-output[[1]]
snplist<-output[[2]]
genes.use<-output[[3]]
vcf_file<-output[[4]]
sets<-output[[5]]
gzs<-output[[6]]
bim<-output[[7]]

genes.use<-intersect(genes.use, gzs)

if(length(genes.use)==0){

     print("All Done!")
     quit(save="no")
}

# Double-check we have only overlapping individuals

indivs.use<-intersect(colnames(rnaseqfile), colnames(vcf_file))

rnaseqfile<-rnaseqfile[which(rnaseqfile[,1] %in% as.character(genes.use)),]
rnaseqfile<-rnaseqfile[,indivs.use]
expdata <- rnaseqfile
expdata<-t(expdata)
dim(expdata)
genes.use<-intersect(genes.use, colnames(expdata))

# get only correct snps in vcf_file

keep<-which(vcf_file$SNP %in% bim[,2])
vcf_file<-vcf_file[keep,]
indivs.use2<-c(colnames(vcf_file)[1:6], indivs.use)
vcf_file<-vcf_file[,indivs.use2]
X<-vcf_file[,-(1:6)]
X <- t(X)

## Now make expdata etc

foldid.use<-X[,1]*0

for(i in 1:10){
    ids<-sets[[i]][,1]
    foldid.use[which(rownames(X)%in%ids)]<-i
}

rm<-c(which(is.na(foldid.use)), which(foldid.use==0))

X<-X[-(rm),]
expdata<-expdata[-(rm),]
foldid.use<-foldid.use[-(rm)]



##### Remove any SNPs which are monomorphic in any corss-validation fold


 uniq.percol<-function(col,X){
         U<-length(unique(X[,col]))
         return(U)
     }


 get_uniques<-function(SET){
    
         Y<-X[-c(sets[[SET]][,1]),]
    
         us<-do.call(c, lapply(1:dim(Y)[2], uniq.percol, Y))
         output<-colnames(X)[which(us==1)]
         return(output)
     }

out<-do.call(c, lapply(1:10, get_uniques))

rm<-unique(out)

if(length(rm)>0){
 X<-X[,-c(which(colnames(X)%in%rm))]
}


#

get.stats<-function(i, cisgenos.full, exppheno){
     return(cv.glmnet(cisgenos.full, exppheno, alpha=0.5, foldid=foldid.use, keep=T))
 }


beta_and_res<-function(GENE, foldid.use){

    genename=gzs[GENE]
    snpinfo<-read.delim(sprintf("/hpc/users/huckil01/snp_lists/%s.onemb.gz", genename), header=F, sep='') # readin snplist
    snpinfo<-snpinfo[which(snpinfo[,1] %in% bim[,2]),]
    cisgenos.full<-X[,which(colnames(X)%in%snpinfo)]
    cisgenos.full <- scale(cisgenos.full, center=T, scale=T)    #scale genotypes
    cisgenos.full[is.na(cisgenos.full)] <- 0

    exppheno<-as.numeric(expdata[,genename])    # rnaseq for gene

    exppheno <- scale(exppheno, center=T, scale=T)  # scale
    exppheno[is.na(exppheno)] <- 0
    best.lam.sim = vector()
    best.cvm.sim = vector()
    pred.matrix = matrix(0,nrow=length(foldid.use),ncol=10)

    glmnet.vec<-vector(10, mode="list")


    glmnet.vec<-mclapply(1:10, get.stats, cisgenos.full, exppheno)  # Run cross-validation step


    for (i in 1:10){

         glmnet.fit<-glmnet.vec[[i]]
         best.lam.ind=which.min(glmnet.fit$cvm)
         best.lam = c(glmnet.fit$cvm[best.lam.ind], glmnet.fit$lambda[best.lam.ind], glmnet.fit$glmnet.fit$df[best.lam.ind], best.lam.ind)
         cvm.best = best.lam[1]
         nrow.max = best.lam[4]
         best.lam.sim[i] = nrow.max
         best.cvm.sim[i] = cvm.best
         pred.matrix[,i] = glmnet.fit$fit.preval[,nrow.max]
     }

    cvm.avg = mean(best.cvm.sim) # average cvm
    nrow.max = as.integer(mean(best.lam.sim)) # best lambda over cv bootstraps

    ret <- as.data.frame(glmnet.fit$glmnet.fit$beta[,nrow.max])
    ret[ret == 0.0] <- NA
    ret.vec = as.vector(ret[which(!is.na(ret)),]) # vector of non-zero betas
    names(ret.vec) = snpinfo[which(!is.na(ret))]
    if(length(ret.vec)==0 && sum(best.lam.sim)==10){ # If we can't make a good enough model, do next-best model

        for (i in 1:10){

            glmnet.fit<-glmnet.vec[[i]]
            best.lam.ind=which.min(glmnet.fit$cvm)+1
            best.lam = c(glmnet.fit$cvm[best.lam.ind], glmnet.fit$lambda[best.lam.ind], glmnet.fit$glmnet.fit$df[best.lam.ind], best.lam.ind)
            cvm.best = best.lam[1]
            nrow.max = best.lam[4]
            best.lam.sim[i] = nrow.max
            best.cvm.sim[i] = cvm.best
            pred.matrix[,i] = glmnet.fit$fit.preval[,nrow.max]
        }

        cvm.avg = mean(best.cvm.sim) # average cvm
        nrow.max = as.integer(mean(best.lam.sim)) # best lambda over cv bootstraps

        ret <- as.data.frame(glmnet.fit$glmnet.fit$beta[,nrow.max])
        ret[ret == 0.0] <- NA
        ret.vec = as.vector(ret[which(!is.na(ret)),]) # vector of non-zero betas
        names(ret.vec) = snpinfo[which(!is.na(ret))]

        write(genename, sprintf("%s/ret.vec.Full.noindel", dir), append=T)

    }

    bestbetas<-ret.vec

    if(length(bestbetas) > 0){
         pred.avg <- rowMeans(pred.matrix)
         res <- summary(lm(exppheno~pred.avg))
         resultsarray<-1:7*0
     #    genename <- as.character(gencode[gene,6])
         resultsarray[1] <- genename
         resultsarray[2] <- cvm.best ###add mean minimum cvm (cross-validated mean-squared error) to results
         resultsarray[3] <- nrow.max ###add mean of best lambda iteration to results
         resultsarray[4] <- glmnet.fit$glmnet.fit$lambda[nrow.max] ###add best lambda to results
         resultsarray[5] <- length(bestbetas) ###add #snps in prediction to results
         resultsarray[6] <- res$r.squared ###lm R2
         resultsarray[7] <- res$coefficients[2,4] ###lm p-value

         ### output bestbetas for PrediXcan
         bestbetalist <- names(bestbetas)
         bestbetainfo <- bim[bestbetalist,]
         betatable<-as.matrix(cbind(bestbetainfo,bestbetas))
         betafile<-cbind(genename, betatable[,2],betatable[,6],betatable[,5],betatable[,7]) ###middle                  colnames(betafile) <- c("Gene","SNP","eff.allele", "ref.allele","beta")
         names(resultsarray)<- c("gene","mean.cvm","mean.lambda.iteration","lambda.min","n.snps","R2","pval")
         write.table(betafile, sprintf("noambig/betas_%s_noambig_1mb", genename), quote=F,row.names=F,sep="\t")
         write.table(resultsarray, sprintf("noambig/res_%s_noambig_1mb", genename), quote=F,row.names=F,sep="\t")
     }

    return(resultsarray)


}

gzs<-intersect(gzs, genes.use)

mclapply(1:length(gzs), beta_and_res, foldid.use)
