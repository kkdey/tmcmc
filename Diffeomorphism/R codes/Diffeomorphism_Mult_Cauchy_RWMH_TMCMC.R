

################   Diffeomorphism Multivariate RWMH and TMCMC  #####################


library(mcmc)
getAnywhere("metrop")
methods(metrop)

#######   Multivariate t  (with diffeomorphism RWMH/TMCMC) ##################


d=50;  ##  dimension of the simulated variable
nsamples=20000;##  sample size in eahc iteration/ replication
L=30; ###   the number of replications we use for finding KS statistic
Mult_General=array(0,c(2,L,nsamples,d));


mu_target=rep(0,d);
Sigma_target=diag(0.7,d)+0.3*rep(1,d)%*%t(rep(1,d));
library(mnormt)
library(fMultivar)


pdf = function(x) 
{
  out=dmt(x,mean=mu_target,Sigma_target,log=TRUE);
  sum=0;
  k=1;
  while(k<=d)
  {
      sum=sum+log(as.numeric(dcauchy2d(x[k],x[k+1],rho=0.3)));
      k=k+2;
      #print(sum)
  }
  return(out)
}


b1 <- morph(b=1);
pdf1<- b1$lud(function(x) pdf(x))

mu1=array(0,d); ###  mean of the proposal density (usually taken to be 0)
propvar=(2.4^2)/d;

# propvar=0.01;
Sigma1=sqrt(propvar)*diag(1,d);


prop_pdf=function(x) dmvnorm(x,mu1,Sigma1);



#####################   RWMH  with diffeomorphism #######################

library(mvtnorm)


for ( l in 1:L)
{
  smpl=matrix(0,nsamples,d);
  smpl[1,]=rnorm(d,10,1);
  i=2;
  counter1=0;
  while (i<=nsamples)
  {
    epsilon= rmvnorm(1,mu1,Sigma1);
    y= smpl[(i-1),]+epsilon;
    accept_rate= min(1,exp(pdf1(y)-pdf1(smpl[(i-1),])));
    w=runif(1);
    if (w<accept_rate)
    {
      smpl[i,]=y;
      i=i+1;
      counter1=counter1+1;
    }
    else  {
      smpl[i,]=smpl[(i-1),];
      i=i+1;   }
  }
  # General[1,l,,]=apply(smpl,c(1,2),b1$inverse);
  Mult_General[1,l,,]=smpl;
  
}

#####################  TMCMC  with diffeomorphism  ######################



for ( l in 1:L)
{
  smpl=matrix(0,nsamples,d);
  smpl[1,]=rnorm(d,10,1);
  i=2;
  counter2=0;
  while (i<=nsamples)
  {
    epsilon= abs(rnorm(1,mu1[1],Sigma1[1,1]));
    y=array(0,d);
    for (j in 1:d)
    {
      w=runif(1)
      if(w<0.5)
      {y[j]= smpl[(i-1),j]+epsilon;}
      else {y[j]= smpl[(i-1),j]-epsilon;}
    }
    
    accept_rate= min(1,exp(pdf1(t(y))-pdf1(t(smpl[(i-1),]))));
    w=runif(1);
    if (w<accept_rate)
    {
      smpl[i,]=y;
      i=i+1;
      counter2=counter2+1;
    }
    else  {
      smpl[i,]=smpl[(i-1),];
      i=i+1;   }
  }
  
  # General[2,l,,]=apply(smpl,c(1,2),b1$inverse);
  Mult_General[2,l,,]=smpl;
  
}


#############   Computing the KS  Statistic  compare ######################


KSval_TMCMC=array(0,nsamples);
KSval_MCMC=array(0,nsamples);

for(d in 1:40)
{
  for(n in 1:nsamples)
  {
    simulate.vec=b1$transform(rcauchy(L,location=0,scale=1));
    KSval_TMCMC[n]=ks.test(Mult_General[2,,n,d],simulate.vec)$statistic;
    KSval_MCMC[n]=ks.test(Mult_General[1,,n,d],simulate.vec)$statistic;
    
  }
  
  
  
  plot(1:nsamples,KSval_TMCMC,col="red",type="l",lwd=1,pch=2,xlab="",ylab="")
  lines(1:nsamples,KSval_MCMC,col="blue",lwd=1,pch=3)
  title(xlab="Time step of run");
  title(ylab="KS test distance");
  title(main="KS plot comparison");
  legend("topright",c("TMCMC","RWMH"),fill=c("red","blue"),border="black");
}