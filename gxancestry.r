library(GxEMM)
##read GCTA grm
source("make_K.R")
#ReadGRMBin=function(prefix, AllN=F, size=4){
# sum_i=function(i){
# return((1+i)*i/2)
# }
#
#  BinFileName=paste(prefix,".grm.bin",sep="")
#  NFileName=paste(prefix,".grm.N.bin",sep="")
#  IDFileName=paste(prefix,".grm.id",sep="")
#  id = read.table(IDFileName)
#  n=dim(id)[1]
#  BinFile=file(BinFileName, "rb");
#  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
#  NFile=file(NFileName, "rb");
#  if(AllN==T){
#    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
#  }
#  else N=readBin(NFile, n=1, what=numeric(0), size=size)
#  i=sapply(1:n, sum_i)
#  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
#}

eur=read.csv("/project2/xuanyao/marie/E-GEUV-1/EUR.QC.TPM.csv",row.names=1,header=T)
yri=read.csv("/project2/xuanyao/marie/E-GEUV-1/Yoruba.QC.TPM.csv",row.names=1,header=T)
all.sample=c(colnames(eur),colnames(yri))
yri.covar=read.table("/project2/xuanyao/xuanyao/trans_ethnic/var_partition/YRI.covar.txt",header=T)
eur.covar=read.table("/project2/xuanyao/xuanyao/trans_ethnic/var_partition/EUR.covar.txt",header=T)
all.covar=rbind(eur.covar,yri.covar)
load("/project2/xuanyao/xuanyao/trans_ethnic/var_partition/eur.exp.Rdata")
load("/project2/xuanyao/xuanyao/trans_ethnic/var_partition/yri.exp.Rdata")
all.ex=rbind(eur.ex,yri.ex)
ldak_loc  <- "./gxemm/ldak5.linux"


dat.yri=read.csv("/project2/xuanyao/xuanyao/trans_ethnic/var_partition/YRI.PIP.BH.csv",header=T)
dat.eur=read.csv("/project2/xuanyao/xuanyao/trans_ethnic/var_partition/EUR.PIP.BH.csv",header=T)
dat=merge(dat.yri,dat.eur,by=8)
output="gxe_he.txt"
write(c("gene","h2_hom","h2_het","h2_original"),file=output,ncol=3,sep="\t")
for(i in 1:6000){
	tryCatch({
	gname=gnames[i]
	cat(i)
	cat("\n")
	cat(paste(gname,"\n",sep=""))
	#gname="ENSG00000185201"
	grm=ReadGRMBin(paste("/project2/xuanyao/marie/E-GEUV-1/LDSC/GCTA/",gname,sep=""))
	K       <- grm$GRM

	fam=read.table(paste("/project2/xuanyao/marie/E-GEUV-1/LDSC/GCTA/",gname,".EandY.fam",sep=""))
	sample.flag=match(fam[,1],all.sample)
	eur.flag=which(is.element(fam[,1],colnames(eur)))


	gene.flag=which(is.element(gnames,gname))

	y=all.ex[sample.flag,gene.flag]
	X=all.covar[sample.flag,]
	Z1=rep(0,nrow(fam))
	Z1[eur.flag]=1
	Z=cbind(Z1,1-Z1)
	X=cbind(X,Z1)
	#out_hom		<- GxEMM( y, X, K, Z, gtype='hom', ldak_loc=ldak_loc)#"mode="HE""
	#out_iid		<- GxEMM( y, X, K, Z, gtype='iid', ldak_loc=ldak_loc) ### need to add etype='iid' for non-discrete environments
	out_hom		<- GxEMM_HE( y, X, K, Z, gtype='hom')#"mode="HE""
	out_iid		<- GxEMM_HE( y, X, K, Z, gtype='iid') ### need to add etype='iid' for non-discrete environments

	### test whether there is any heritability assuming the Hom model
	#Waldtest( out_hom$h2, out_hom$h2Covmat[1,1] )   

	### test for genetic heterogeneity using IID model, which assumes that h2 is equal across all environments
	#Waldtest( out_iid$h2[2], out_iid$h2Covmat[2,2] )

	### tests for genetic heterogeneity using Free model
	#pvalue=MVWaldtest( out_free$sig2s[2:3], out_free$sig2Var[2:3,2:3] ) 
	filename=paste("./results_he/",gname,".Rdata",sep="")
	write(c(gname,out_iid$h2[1],out_iid$h2[2],out_hom$h2[1]),file=output,ncol=3,sep="\t",append=TRUE)
	out_free	<- GxEMM_HE( y, X, K, Z, gtype='free', etype='free')
	save( out_hom, out_iid, file=filename)

},error=function(e){})
}
