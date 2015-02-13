function Out = U(q, mu, Sigma)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
pdf=@(x)mvnpdf(x,mu,Sigma); 
Out = -log(pdf(q));
end

