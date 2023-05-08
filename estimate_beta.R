estimate_beta= function(BETA,SE,R,lx=100,est=TRUE) {

weightx = matrix(NA, dim(BETA)[1], dim(BETA)[2])

if(est==TRUE){
for (SS in 1)
{
S = diag(SE[,SS])
SRS = (S%*%R%*%S)
eig = eigen(SRS)
d = eig$values
U = eig$vectors
indexd = which(d < 0)
if (length(indexd) > 0) d[indexd] = 1
d2 = 1/sqrt(d)
if (length(indexd) > 0) d2[indexd] = 0
SRS12 = U%*%diag(d2)%*%t(U)
yg = c(SRS12%*%BETA[,SS])
Mg = as.matrix(SRS12%*%(S%*%R%*%diag(1/SE[,SS])))
fitw = cv.glmnet(Mg,yg,intercept=FALSE,standardize=FALSE,alpha=0,nlambda=lx,family="gaussian")
weightx[,SS] = c(coef(fitw,s="lambda.min")[-1,1])
}
}else{
weightx=as.matrix(BETA[,1])
}
return(weightx)
}