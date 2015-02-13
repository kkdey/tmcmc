function Out= MCMC_fun(nsamples, d, x0, pdf, Sigma1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
i=2;
counter=0;
smpl=zeros(nsamples,d);
smpl(1,:)=x0;

while i<=nsamples
    epsilon= mvnrnd(zeros(1,d), Sigma1,1);
    y= smpl((i-1),:)+epsilon;
    if pdf(y)~=0
        accept_rate= min(1,(pdf(y)/pdf(smpl((i-1),:))));
        w=rand(1);
        if w<accept_rate
            smpl(i,:)=y;
            i=i+1;
            counter=counter+1;
        else
            smpl(i,:)=smpl((i-1),:);
            i=i+1;
            
        end
    end
    
    
end

Out=smpl;

end

