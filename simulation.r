####simulation for rg=1, demonstraing biases due to constraints#####
####real genotype, wide range of heritabiity and rg=1
####Use GCTA to estimate rg ######
library(MASS)
require("snpStats")
set.seed(1234)

standardize=function(geno){

	snp_mean=apply(geno,2,mean,na.rm=T)
	for(i in 1:nrow(geno)){
		geno[i,]=geno[i,]-snp_mean
	}
	snp_sd=apply(geno,2,sd,na.rm=T)
	for(i in 1:ncol(geno)){
		geno[,i]=geno[,i]/snp_sd[i]
	}
	geno[which(is.na(geno))]=0
	return(geno)

}


simulate_betaX=function(geno1,geno2,p,h2,rg){
        M=ncol(geno1)
        N1=nrow(geno1)
	N2=nrow(geno2)
	s_causal=rbinom(M,1,p)
	sigma=matrix(c(1,rg,rg,1),ncol=2,nrow=2)
        beta=mvrnorm(sum(s_causal!=0),rep(0,2),sigma)*sqrt(h2/(M*p))
	sflag=which(s_causal==1)
        G1=geno1[,sflag]%*%beta[,1]
        G2=geno2[,sflag]%*%beta[,2]
        noise1=rnorm(N1,0,sqrt(1-h2))
        noise2=rnorm(N2,0,sqrt(1-h2))
	phe1=G1+noise1
	phe2=G2+noise2
        E=c(phe1,phe2)
        return(E)
}


argv=commandArgs(trailingOnly=T)
#input h2 and rg
h2=as.numeric(argv[1])
rg=as.numeric(argv[2])

p_causal=as.numeric(argv[3])
input_dir=argv[4]
output_dir=argv[5]
subdir=paste0("h2_",h2)
constraint_dir="./constraint/"

dir.create(file.path(input_dir,subdir),recursive=TRUE)
dir.create(file.path(output_dir,subdir),recursive=TRUE)
dir.create(file.path(constraint_dir,subdir),recursive=TRUE)
####load genotype######
gene="ENSG00000129116"
YRI=read.table("../var_partition/YRI.sampleID.txt")
snps <- read.plink("/project2/xuanyao/marie/E-GEUV-1/LDSC/GCTA/ENSG00000129116.EandY")
geno=as(snps$genotypes,"numeric")
sample=snps$fam
flag=which(is.element(sample$member,YRI[,1]))
N1=length(flag)
N2=nrow(geno)-N1
geno1=geno[-flag,] ##EUR genotypes
geno2=geno[flag,] ##YRI genotypes

geno1=standardize(geno1)
geno2=standardize(geno2)
sample1=sample[-flag,]
sample2=sample[flag,]
sample_new=rbind(sample1,sample2)
for(i in 1:10000){
	exp=simulate_betaX(geno1,geno2,p_causal,h2,rg)
	exp2=cbind(c(exp[1:N1],rep("NA",N2)),c(rep("NA",N1),exp[-c(1:N1)]))
	write.table(cbind(sample_new[,1:2],exp2),paste0(input_dir,subdir,"/pheno_iter",i,"_h2_",h2,"_rg_",rg,"_p_",p_causal,".txt"),row.names=F,col.names=F,quote=F)

	####
	system(paste("gcta64 --reml-bivar --reml-bivar-no-constrain --reml-no-constrain --reml-maxit 100 --grm  /project2/xuanyao/marie/E-GEUV-1/LDSC/GCTA/",gene,"  --pheno ",input_dir,subdir,"/pheno_iter",i,"_h2_",h2,"_rg_",rg,"_p_",p_causal,".txt  --out ",output_dir,subdir,"/iter",i,"_h2_",h2,"_rg_",rg,"_p_",p_causal,"_output",sep=""))
	system(paste("gcta64 --reml-bivar  --reml-maxit 100 --grm  /project2/xuanyao/marie/E-GEUV-1/LDSC/GCTA/",gene,"  --pheno ",input_dir,subdir,"/pheno_iter",i,"_h2_",h2,"_rg_",rg,"_p_",p_causal,".txt  --out ",constraint_dir,subdir,"/iter",i,"_h2_",h2,"_rg_",rg,"_p_",p_causal,"_output",sep=""))

	}
ddir=paste0(input_dir,subdir)
unlink(ddir,recursive=TRUE)
