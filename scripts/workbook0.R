#helper function to save plots 
savePlot <- function(plt, outfile) {
  pdf(outfile)
  print(plt)
  dev.off()
}
savePlot2 <- function(plt,subplt, outfile) {
  pdf(outfile)
  print(plt)
  print(subplt)
  dev.off()
}
###load the data
filename='../data0/sdcate12T.lt'
filename='../data0/sdcate9T.lt'
#filename='../data0/sdcate8T.lt'#bad force stability?
filename='../data0/sdcate6T.lt'
#filename='../data0/sdcate5BT.lt'
#filename='../data0/sdcate5AT_2.lt'
D<-read.table(filename,header=F)

#preprocess the data slightely 
colnames(D)<-c("time","x","force")
N<-length(D$time)#number of points
max_lag<-floor(N/20.0)
lags<-seq(5,max_lag,5)

mean.dX<-array()
var.dX<-array()
tau<-array()
weights<-array()
i<-1
#calculate statistics of displacements:
for(lag in lags){
  sub_idx<-seq(1,N,lag)#subsample
  X=D$x[sub_idx]
  ts=D$time[sub_idx]
  mean.dX[i]=mean(diff(X))
  var.dX[i]=var(diff(X))
  tau[i]=mean(diff(ts))
  weights[i]=length(sub_idx)^0.5
  i=i+1
}

RV.dX=var.dX/(mean.dX^2/tau^2)
weights=weights/sum(weights)

stats=data.frame(tau,mean.dX,var.dX,RV.dX)

##basic fits:
source("./step_stats.R")
#
fit.mean<-lm(mean.dX ~ tau -1,data=stats)
summary(fit.mean)
#fit variance
sigma0<-sqrt(0.5*var.dX[1])##initial guess for noise 
v0<-as.numeric(coef(fit.mean)["tau"])
A0<-mean(var.dX[1:5])/v0^2;
k0<-1.0
B0<-1.0
#NLS fit:
fit.rv<-nls( RV.dX ~ var.model( tau, A, B, k ), 
              data=stats, 
              start=list(A=A0,B=B0,k=k0),
              #lower=list(A=0.000,B=0,k=1e-6),
              weights=weights,
              control=nls.control(maxiter = 1200, tol = 1e-04, minFactor = 1e-6,
                                  printEval = FALSE, warnOnly = FALSE));
summary(fit.rv)

cat(filename)
cat('\n')
cat(sprintf('size:=%d, max_lag:=%d',N,max_lag))
cat('\n')

k<-summary(fit.rv)$coefficients[3,1]
B<-summary(fit.rv)$coefficients[2,1]
B.stderr<-summary(fit.rv)$coefficients[2,2]
A<-summary(fit.rv)$coefficients[1,1]
pmo<-B/2.0
pmo.stderr<-B.stderr/2.0
nu<-1.0/(1.0+pmo) ##prob. of the moving state
nu.stderr<-0.5*nu*(B.stderr/B)
Vmax<-as.numeric(fit.mean$coeff['tau'])/nu
V<-as.numeric(fit.mean$coeff['tau'])
V.stderr <- sqrt(diag(vcov(fit.mean)))[1]

##parameter covariance
cov<-vcov(fit.rv)
##noise amplitude
sigma_noise<-sqrt(0.5*abs(A)*V^2)

cat(sprintf("PMO: %f %f\n",pmo,pmo.stderr))
cat(sprintf("PO: %f  %f\n",1.0-nu,nu.stderr/nu*(1.0-nu)))
cat(sprintf("V:  %f  %f\n",V,V.stderr[1]))
cat(sprintf("Vmax:  %f %f\n",Vmax, Vmax*nu.stderr))

cat(sprintf("pause-length(sec)=%2.4f\n",(1.0+pmo)/k))
cat(sprintf("pause-freq(1/sec):=%2.4f\n",nu*k*pmo/(1.0+pmo)))

#basic plots
library(ggplot2)
library(grid)
#plot var:
RV.dX.predict<-data.frame(tau,RV.dX=predict(fit.rv,tau))
#
var.plot<-qplot(tau,RV.dX,data=stats,xlab="tau(sec)",ylab="RV(tau)",main=sprintf("file:=%s,rescaled variance",filename))
var.plot<-var.plot+geom_point(size=2.5, colour="#CC0000")
var.plot<-var.plot+geom_line(data=RV.dX.predict,aes(tau,RV.dX),linetype=2)
print(var.plot)

#
vp<-viewport(width = 0.40, height = 0.35, x = 0.30, y = 0.7)
trace.plot<-ggplot(D,aes(time,x),)+geom_line(color = "black")+ylab("x(nt)")+xlab("time(sec)")+ggtitle("position")
print(trace.plot,vp=vp)

save_plot=T

if(save_plot){
 pdf_file=sprintf("%s.pdf",filename)
 pdf(pdf_file)
 print(var.plot)
 print(trace.plot,vp=vp)
 dev.off()
 pdf_file=sprintf("%s_hist_force.pdf",filename)
 pdf(pdf_file)
 force.plot<-ggplot(D,aes(x=force))+geom_histogram(colour="green",fill="blue")
 print(force.plot)
 dev.off()
}



