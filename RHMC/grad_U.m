function Out = grad_U(q,mu,Sigma)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Out=inv(Sigma)*(q-mu)';
end

