generate_fold_ids <- function(n_samples, n_folds=10) {
  n <- ceiling(n_samples / n_folds)
  fold_ids <- rep(1:n_folds, n)
  sample(fold_ids[1:n_samples])
}


library(glmnet)
set.seed(123)

eur=read.csv("/project2/xuanyao/marie/E-GEUV-1/EUR.QC.TPM.csv",row.names=1,header=T)
yri=read.csv("/project2/xuanyao/marie/E-GEUV-1/Yoruba.QC.TPM.csv",row.names=1,header=T)
all.tpm=merge(eur,yri,by=0)
library(preprocessCore)
library(phenix)

eur=all.tpm[,2:359]
yri=all.tpm[,360:446]
all.eur <- normalize.quantiles(as.matrix(eur))
all.eur.t=t(all.eur)
all2.eur=quantnorm(all.eur.t)
all.eur=t(all2.eur)
ex.eur=t(all.eur)

all.yri <- normalize.quantiles(as.matrix(yri))
all.yri.t=t(all.yri)
all2.yri=quantnorm(all.yri.t)
all.yri=t(all2.yri)
ex.yri=t(all.yri)

gnames=as.character(all.tpm[,1])

n_folds=10
covar.yri=read.table("../var_partition/YRI.covar.txt",header=T)
covar.eur=read.table("../var_partition/EUR.covar.txt",header=T)
#write(c("gene","r2.yri","r2.self","i"),file="elasticnet_pred.txt",ncol=4,sep="\t")
for(i in 6865:length(gnames)){
# Fit on all dataA
	gene=gnames[i]
	if(file.exists(paste("/project2/xuanyao/marie/E-GEUV-1/finemap/varLD/",gene,".varLDgenotype.EUR.csv",sep=""))){
	eur_geno_file=read.table(paste("/project2/xuanyao/marie/E-GEUV-1/finemap/varLD/",gene,".varLDgenotype.EUR.csv",sep=""))	
	yri_geno_file=read.table(paste("/project2/xuanyao/marie/E-GEUV-1/finemap/varLD/",gene,".varLDgenotype.Yoruba.csv",sep=""))	
        expression.yri=ex.yri[,i]
	reg.yri=lm(expression.yri~as.matrix(covar.yri))
	adj_yri=reg.yri$residuals[-1]
	expression.eur=ex.eur[,i]
	reg.eur=lm(expression.eur~as.matrix(covar.eur))
	adj_eur=reg.eur$residuals[-1]

	eur_geno=t(as.matrix(eur_geno_file[,-c(1,2)])-1)
	flag=sort(sample(nrow(eur_geno),271))
	eur_geno_sub=eur_geno[flag,]
	adj_eur=adj_eur[flag]

	yri_geno=t(as.matrix(yri_geno_file[,-c(1,2)])-1)
	cv_fold_ids <- generate_fold_ids(length(adj_eur), n_folds)
	fit <- tryCatch(cv.glmnet(eur_geno_sub, adj_eur, nfolds = n_folds, alpha = 0.5, type.measure='mse', foldid = cv_fold_ids, keep = TRUE),
                      error = function(cond) {message('Error'); message(geterrmessage()); list()})

	best_lam_ind <- which.min(fit$cvm)
        if (fit$nzero[best_lam_ind] > 0) {
	cat(i)
	cat("\n")

         	weights <- fit$glmnet.fit$beta[which(fit$glmnet.fit$beta[,best_lam_ind] != 0), best_lam_ind]
		pred=as.matrix(yri_geno[,which(fit$glmnet.fit$beta[,best_lam_ind] != 0)])%*%weights
		pred.self=as.matrix(eur_geno[-flag,which(fit$glmnet.fit$beta[,best_lam_ind] != 0)])%*%weights
		
		r2.self=(cor(expression.eur[-c(1,flag+1)],pred.self))^2
		r2=(cor(expression.yri[-1],pred))^2
		summ=c(gene,r2,r2.self,i)
		write(summ,file="elasticnet_pred.txt",append=TRUE,ncol=4,sep="\t")
		}
	}#end if

}
