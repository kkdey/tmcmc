function [out] = Select_temp(minbeta, mu1, mu2, Sigma11, Sigma12, p)
% Selecting the inverse temperatures for swapping 
%   Stochastic approximation
L=20; % length of each iteration
current_beta=1;
beta_array=1;
while current_beta > minbeta
    rho=zeros(1,L+1);
    newbeta=zeros(1,L);
    alpha = zeros(1,L);
    for n=1:L
        newbeta(n)= current_beta*inv(1+exp(rho(n)));
        [x_curr1 x_curr2] = Rand_generate(current_beta, mu1, mu2, Sigma11, Sigma12, p);
        [x_new1 x_new2] = Rand_generate(newbeta(n), mu1, mu2, Sigma11, Sigma12, p);
        x_curr=[x_curr1 x_curr2];
        x_new= [x_new1 x_new2];
        hand=@(x,y) pdfmix([x y],mu1,mu2,Sigma11,Sigma12,p);
        hand_log = @(x,y) log(hand(x,y));
        b= -(newbeta(n)-current_beta)*(hand_log(x_new1,x_new2)-hand_log(x_curr1,x_curr2));
        alpha(n)= min(1,b);
        rho(n+1)= rho(n)+ (1/n)*(alpha(n)-0.234);
    end
    
    current_beta= current_beta*inv(1+exp(rho(L+1)));
    beta_array=[beta_array current_beta];
end

out=beta_array;
end

