function [out1, out2] = Rand_generate(beta,mu1,mu2,Sigma11,Sigma12,p)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
hand=@(x,y) pdfmix([x y],mu1,mu2,Sigma11,Sigma12,p);
beta_hand= @(x,y) power(hand(x,y),beta);
% fplot(hand,[-2 12]);
d=length(mu1);
x0=zeros(1,d);
nsamples=4000;
smpl_TMCMC=zeros(nsamples,d);
sigmaprop=power(10,2)/d;
smpl_TMCMC(1,:)=x0;
i=2;
counterT=zeros(1,nsamples);
counter=0;

while i<=nsamples
    epsilon= abs(normrnd(0,sqrt(sigmaprop)));
    y=zeros(1,d);
    for x=1:d
        u=rand(1);
        if u<0.5
            y(x)=smpl_TMCMC((i-1),x)+epsilon;
        else
            y(x)=smpl_TMCMC((i-1),x)-epsilon;
        end
    end
    
    accept_rate= min(1,(beta_hand(y(1),y(2))/beta_hand(smpl_TMCMC((i-1),1),smpl_TMCMC((i-1),2))));
    w=rand(1);
    if w<accept_rate
        smpl_TMCMC(i,:)=y;
        counterT(i)=1;
        counter=counter+1;
        i=i+1;
        
    else
        smpl_TMCMC(i,:)=smpl_TMCMC((i-1),:);
        counterT(i)=0;
        i=i+1;
       
    end
end
l=10;
r=smpl_TMCMC((nsamples-l),:);
while (r(1)<-3)||(r(1)>13)
    while (r(2)<-3) || (r(2)>13)
        l=l+1;
        r=smpl_TMCMC((nsamples-l),:);
    end
    
end

out1=r(1);
out2=r(2);



end

