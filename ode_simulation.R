closed.sir.model <- function (t, x, params) {
  ## first extract the state variables
  S <- x[1]
  L1 <-x[2]
  L2 <-x[3]
  I <- x[4]
  R <- x[5]
  ## now extract the parameters
  pi <-params["pi"]
  beta <- params["beta"]
  mu <- params["mu"]
  mut <- params["mut"]
  p <- params["p"]
  x <- params["x"]
  f <- params["f"]
  q <- params["q"]
  r <- params["r"]
  y <- params["y"]
  k <- params["k"]
  a <- params["a"]
  z <- params["z"]
  N <- params["N"]
  
  ## now code the model equations
  dSdt <- pi*N-beta*I*S/N-mu*S+r*R
  dL1dt <- (1-p)*beta*I*S/N-(f+ x*beta*I/N + mu)*L1
  dL2dt <- (1-q)*x*beta*I*L1/N-(y+mu)*L2+ z*R
  dIdt <- p*beta*I*S/N + (f+x*q*beta*I/N)*L1 + y*L2 -(k+mu+mut)*I+a*R
  dRdt <- k*I-(r+mu+a+z)*R 
 ## combine results into a single vector
  dxdt <- c(dSdt,dL1dt,dL2dt,dIdt,dRdt) 
  ## return result as a list!
  list(dxdt)                              
}

params <- c(pi=0.02,beta=30,mu=0.02,mut=0.075,p=0.005,x=0.2,z=0.97,f=0.006,q=0.1,r=0.0,y=0.0015,k=0.68,a=0.03,N=100001)

times <- seq(from=0,to=300,by=0.5) # returns a sequence
xstart <- c(S=100000,L1=0.00,L2=0.00,I=1.00,R=0.00)     # initial conditions

require(deSolve)

out <- as.data.frame(
                     ode(
                         func=closed.sir.model,
                         y=xstart,
                         times=times,
                         
                         parms=params
                         )
                     )
png("/home/chacha/Downloads/ecr-figuresR/TL30.png")
plot(S~time,data=out,type='l',col="green",lwd=4,xlab="Time (years)",ylab="Population",ylim=c(0,100000),cex.lab=1.3)
lines(L1~time,data=out,type='l',col="blue",lwd=4)
lines(L2~time,data=out,type='l',col="brown",lwd=4)
lines(I~time,data=out,type='l',col="red",lwd=4)
lines(R~time,data=out,type='l',col="black",lwd=4)

legend("topright",legend=c("S","L1","L2","I","T"),col=c("green","blue","brown","red","black"),lty=1,cex=0.85,lwd=4)
dev.off()
