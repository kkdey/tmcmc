function [out] = TMC3( mu1, mu2, Sigma11, Sigma12,p )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
minbeta=0.2;
beta_set=Select_temp(minbeta,mu1,mu2,Sigma11,Sigma12,p);
d_beta=length(beta_set);
TMC3=cell(1,d_beta);
d=length(mu1);

x0=ones(1,d);

for k=1:d_beta
    TMC3{1,k}(1,:)=x0;
end
nsamples=4000;
i=2;
propvar=power(2.4,2)/d;
sigmaprop=propvar;
Sigma1=propvar*eye(d);
mu1=zeros(1,d);

counter=0;

while i<=nsamples
    %smpl=zeros(nsamples,d);
    for k=1:d_beta
        epsilon= abs(normrnd(mu1(1),sqrt(sigmaprop)));
        y=zeros(1,d);
        for x=1:d
            u=rand(1);
            if u<0.5
                y(x)=TMC3{1,k}((i-1),x)+epsilon;
            else
                y(x)=TMC3{1,k}((i-1),x)-epsilon;
            end
        end
        
        hand=@(x) pdfmix(x,mu1,mu2,Sigma11,Sigma12,p);
        beta_hand= @(x) power(hand(x),k);
        accept_rate= min(1,(beta_hand(y)/beta_hand(TMC3{1,k}((i-1),:))));
        w=rand(1);
        if w<accept_rate
            TMC3{1,k}(i,:)=y;
            counter=counter+1;
        else
            TMC3{1,k}(i,:)=TMC3{1,k}((i-1),:);
        end
    end
    
    for k=2:d_beta
        beta_hand_1= @(x) power(hand(x),k);
        beta_hand_2= @(x) power(hand(x),k-1);
        u= TMC3{1,k}(i,:);
        v= TMC3{1,k-1}(i,:);
        swap_rate= min(1,(beta_hand_1(v)*beta_hand_2(u)/beta_hand_1(u)*beta_hand_2(v)));
        w=rand(1);
        if w<swap_rate
            temp=TMC3{1,k}(i,:);
            TMC3{1,k}(i,:)=TMC3{1,k-1}(i,:);
            TMC3{1,k-1}(i,:)=temp;
        end
    end
    
    i=i+1;
end

out=TMC3{1,1};

    
end


