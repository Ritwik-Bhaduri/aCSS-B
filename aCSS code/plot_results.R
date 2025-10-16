#If you run run_sim.R 2000 times, for run = 1,...,500 and example_num = 1,...,4, this script will produce the exact plots from the paper from the saved outputs.

dir = getwd() # Set this to wherever you saved the results, i.e., save.dir in RC_final.R
width = 4
height = 4
alpha=0.05
repeats = 500

#########################################################################
# Multivariate t example

deltas = seq(from=0,to=8,by=2)
pwrs = rep(0,length(deltas))
pwrs.o = rep(0,length(deltas))
pwrs.s = rep(0,length(deltas))
for(i in 1:repeats){
	filename = paste0(dir,"/MT/MTr",i,".rdata")
	load(filename)
	pwrs = pwrs + c(results[,2]<=alpha)/repeats
	pwrs.o = pwrs.o + c(results[,1]<=alpha)/repeats
	pwrs.s = pwrs.s + c(results[,3]<=alpha)/repeats
}
pdf(file=paste0(dir,"/Multivariatet.pdf"),width=width,height=height)
plot(deltas,rep(0,length(deltas)),type="l",col="white",ylab="Power",ylim=c(0,1),xlab="True d.f. - Null d.f.",main="Multivariate t")
points(deltas,pwrs,type="l",col="black",lty=1)
points(deltas,pwrs.o,type="l",col="black",lty=2)
points(deltas,pwrs.s,type="l",col="black",lty=3)
abline(h=0.05,lty=1,col="grey")
legend("topleft",c("aCSS","oracle","score"),col=c("black","black","black"),lty=c(1,2,3))
dev.off()

#########################################################################
# Logistic Regression example

deltas = seq(from=0,to=1,by=0.1)
pwrs = rep(0,length(deltas))
pwrs.o = rep(0,length(deltas))
for(i in 1:repeats){
	filename = paste0(dir,"LR/LRr",i,".rdata")
	load(filename)
	pwrs = pwrs + c(results[,2]<=alpha)/repeats
	pwrs.o = pwrs.o + c(results[,1]<=alpha)/repeats
}
pdf(file=paste0(dir,"/LogisticRegression.pdf"),width=width,height=height)
plot(deltas,rep(0,length(deltas)),type="l",col="white",ylab="Power",ylim=c(0,1),xlab="Coefficient on X",main="Logistic Regression")
points(deltas,pwrs,type="l",col="black",lty=1)
points(deltas,pwrs.o,type="l",col="black",lty=2)
abline(h=0.05,lty=3)
legend("topleft",c("aCSS","oracle"),col=c("black","black"),lty=c(1,2))
dev.off()

#########################################################################
# Behrens-Fisher example

deltas = seq(from=0,to=1,by=0.1)
pwrs = rep(0,length(deltas))
pwrs.o = rep(0,length(deltas))
pwrs.s = rep(0,length(deltas))
for(i in 1:repeats){
	filename = paste0(dir,"BF/BFr",i,".rdata")
	load(filename)
	pwrs = pwrs + c(results[,2]<=alpha)/repeats
	pwrs.o = pwrs.o + c(results[,1]<=alpha)/repeats
	pwrs.s = pwrs.s + c(results[,3]<=alpha)/repeats
}
pdf(file=paste0(dir,"/BehrensFisher.pdf"),width=width,height=height)
plot(deltas,rep(0,length(deltas)),type="l",col="white",ylab="Power",ylim=c(0,1),xlab=expression(mu^(1)-mu^(0)),main="Behrens-Fisher")
points(deltas,pwrs,type="l",col="black",lty=1)
points(deltas,pwrs.o,type="l",col="black",lty=2)
points(deltas,pwrs.s,type="l",col="black",lty=3)
abline(h=0.05,lty=1,col="grey")
legend("topleft",c("aCSS","oracle","score"),col=c("black","black","black"),lty=c(1,2,3))
dev.off()

#########################################################################
# Gaussian Spatial example

deltas = seq(from=0,to=1,by=0.2)
pwrs = rep(0,length(deltas))
pwrs.o = rep(0,length(deltas))
for(i in 1:repeats){
	filename = paste0(dir,"GS/GSr",i,".rdata")
	load(filename)
	pwrs = pwrs + c(results[,2]<=alpha)/repeats
	pwrs.o = pwrs.o + c(results[,1]<=alpha)/repeats
}
pdf(file=paste0(dir,"/GaussianSpatial.pdf"),width=width,height=height)
plot(deltas,rep(0,length(deltas)),type="l",col="white",ylab="Power",ylim=c(0,1),xlab="Anisotropy Parameter",main="Gaussian Spatial")
points(deltas,pwrs,type="l",col="black",lty=1)
points(deltas,pwrs.o,type="l",col="black",lty=2)
abline(h=0.05,lty=3)
legend("topleft",c("aCSS","oracle"),col=c("black","black"),lty=c(1,2))
dev.off()