function [out] = pdfcauchy(x,mu,sigma)
z=(x-mu)/sigma;
pdf1= @(x) 1*(1/pi)*(1/1+x^2);
p=1;
for l=1:length(z)
    p=p*pdf1(z(l));
end
out=p;    
end

%addpath('C:\Documents and Settings\Admin\My Documents\MATLAB\RDMH');