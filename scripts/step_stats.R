#
# algorithm for scoring step intensity in molecular motor trajectory
# send comments to sergey.plyasunov@gmail.com
# Physics Department, UC Berkeley
# first version March 2006
#current version July 31, 2012
#
##basic-stats
kurtosis <- function(x) {
                m4 <- mean((x-mean(x))^4)
                kurt <- m4/(sd(x)^4)-3.0e0
                return(kurt)
}

skewness <-  function(x) {
                m3 <- mean((x-mean(x))^3)
                skew <- m3/(sd(x)^3)
                return(skew)
}

JB<-function(x){

                N<-length(x)
                sk<-skewness(x);
                ku<-kurtosis(x);
                val<-(N/6.0)*( sk^2 + 0.25 *(ku)^2);
                return(val);
}


LaplaceTransform<-function(u,X){

                N<-length(X);
                val<-0.0e0;

                for(i in c(1:N)){
                      #val <- val+ (1.0/N)* cos(u*X[i]);
		      val <- val+ (1.0/N)* sin(u*X[i]);
                 }

                return(val);
}

LaplaceTransform.model<-function( u, v0, p, sigma0 ){

                #v0<-params[1];
                #p<-params[2];
                #sigma0<-params[3];

                #val<-( (1.0-p)*cos(u*v0) + p )*exp(-0.5*u^2*sigma0^2);
		val<-(1.0-p)*sin(u*v0)*exp(-0.5*u^2*sigma0^2);

                return(val);
}
##model for rescaled variance as a function of time-window "tau" and model params:
var.model<-function(tau, A, B,  kappa ){
			
	val = A + B * ( kappa*tau-(1.0 - exp(-kappa*tau)) )/kappa^2;
	return(val)
}

var.model2<-function(tau, A, B , k){
	val<-(A/tau^2 + B *(k*tau+1+exp(-k*tau))/tau^2 )
	return(val)
}
