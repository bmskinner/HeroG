# Use RefFreeEWAS to account for cell composition in WGS data
# Based on tutorial at https://akhilesh362.wordpress.com/

install.packages("RefFreeEWAS")
library(RefFreeEWAS)
# Replace below "..." with the directory containing the data
y = read.table("...\input_data.txt",sep="\t",header=T)
y1 = y[,-1];rownames(y1) = y[,1]
miss = which(is.na(y1),arr.ind=TRUE)[,1]
y1 = y1[-miss,]
edata = as.matrix(log2(y1/(1-y1)))

covar1 = read.table("input_phenotype.txt",header=FALSE);
covar2<-read.table("covariates.txt",header=FALSE)
x1 = as.numeric(covar1[,2])
x2<-as.numeric(covar2[,3])
x3<-as.numeric(covar2[,4])
x4<-as.numeric(covar2[,5])
x5<-as.numeric(covar2[,6])
x6<-as.numeric(covar2[,7])
tmp1 = lm(t(edata)~x1);
dim3 = EstDimRMT(cbind(t(coef(tmp1)), t(resid(tmp1))), FALSE)$dim
test1 = RefFreeEwasModel(y2,cbind(1,x1,x2,x3,x4,x5,x6),dim3)
testBoot1 = BootRefFreeEwasModel(test1, 50)

res1 = summary(testBoot1); est1 = res1[, 2, 1, 'mean']
est.sd1 = res1[, 2, 1, 'sd']
Q21 = (est1/est.sd1)^2; numCov = 1
P01 = 2*(pt(-sqrt(Q21),(204-(dim3+6))))
fdr.pvalue =  p.adjust(P01, method = "fdr", n = length(P01))
i = which(fdr.pvalue <= 0.05)
est_ewas1 = est1[i]; sd_ewas1 = est.sd1[i]; P0_ewas1 = P01[i]
fdr_ewas1 = fdr.pvalue1[i]
cpg_ewas1 = cbind(est_ewas1,sd_ewas1,P0_ewas1,fdr_ewas1)
write.table(cpg_ewas1,file="...\filname.csv",sep=",")