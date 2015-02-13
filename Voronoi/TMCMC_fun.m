function Out = TMCMC_fun( nsamples, d, x0, pdf, sigmaprop )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
smpl_TMCMC=zeros(nsamples,d);
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
    
    accept_rate= min(1,(pdf(y)/pdf(smpl_TMCMC((i-1),:))));
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

Out=smpl_TMCMC;

end

