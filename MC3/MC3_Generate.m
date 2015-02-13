function [out] = MC3_Generate(mu1, mu2, Sigma11, Sigma12, p)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
minbeta=0.2;
beta_set=Select_temp(minbeta,mu1,mu2,Sigma11,Sigma12,p);
d_beta=length(beta_set);
MC3=cell(1,d_beta);
d=length(mu1);

x0=zeros(1,d);

for k=1:d_beta
    MC3{1,k}(1,:)=x0;
end
nsamples=4000;
i=2;
propvar=power(2.4,2)/d;
Sigma1=propvar*eye(d);
mu1=zeros(1,d);

counter=0;

while i<=nsamples
    %smpl=zeros(nsamples,d);
    for k=1:d_beta
        epsilon= mvnrnd(mu1, Sigma1,1);
        y= MC3{1,k}((i-1),:)+epsilon;
        hand=@(x,u) pdfmix([x,u],mu1,mu2,Sigma11,Sigma12,p);
        beta_hand= @(x,u) power(hand(x,u),k);
        accept_rate= min(1,(beta_hand(y(1),y(2))/beta_hand(MC3{1,k}((i-1),1),MC3{1,k}((i-1),2))));
        w=rand(1);
        if w<accept_rate
            MC3{1,k}(i,:)=y;
            counter=counter+1;
        else
            MC3{1,k}(i,:)=MC3{1,k}((i-1),:);
        end
    end
    
    for k=2:d_beta
        beta_hand_1= @(x,y) power(hand(x,y),k);
        beta_hand_2= @(x,y) power(hand(x,y),k-1);
        u= MC3{1,k}(i,:);
        v= MC3{1,k-1}(i,:);
        swap_rate= min(1,(beta_hand_1(v(1),v(2))*beta_hand_2(u(1),u(2))/beta_hand_1(u(1),u(2))*beta_hand_2(v(1),v(2))));
        w=rand(1);
        if w<swap_rate
            temp=MC3{1,k}(i,:);
            MC3{1,k}(i,:)=MC3{1,k-1}(i,:);
            MC3{1,k-1}(i,:)=temp;
        end
    end
    
    i=i+1;
end

out=MC3{1,1};

    
end

