set.seed(1)
nrep = 10000
n = 100
M = 500
rho = 0.97

pval_parboot = pval_css = rep(0,nrep)
for(rep in 1:nrep){
	X = rnorm(n)
	Y = rnorm(n)
	Z = Y * rho + rnorm(n) * sqrt(1-rho^2)
	thetah = sum(Z*X)/sum(Z^2)
	Xt_parboot = matrix(rnorm(n*M),n,M) + thetah * Z
	Xt_css = (diag(n)-Z%*%t(Z) / sum(Z^2))%*%matrix(rnorm(n*M),n,M) + thetah * Z

	TX =  sum(X*Y)^2 / (sum(X*Z))^2
	TXt_parboot =  (t(Xt_parboot)%*%Y)^2 / (t(Xt_parboot)%*%Z)^2
	TXt_css = (t(Xt_css)%*%Y)^2 / (t(Xt_css)%*%Z)^2

	pval_parboot[rep] = (1+sum(TXt_parboot>=TX))/(1+M)
	pval_css[rep] = (1+sum(TXt_css>=TX))/(1+M)
}

par(mfrow=c(1,2))
hist(pval_parboot,xlim=c(0,1),main='Parametric bootstrap',xlab='p-value')
hist(pval_css,xlim=c(0,1),main='Co-sufficient sampling',xlab='p-value')
