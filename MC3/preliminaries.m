% MCMCMC method and its TMCMC analog

fplot(@(x) pdfmix(x,0,10,1,1,1),[-4 4])
xlabel('x');
ylabel('density');

handle=@(x) pdfmix(x,0,8,1,1,0.5);
beta=0;
fplot(@(x) power(handle(x),beta),[-4 12]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu1=repmat(0,1,d);
mu2=repmat(4,1,d);
Sigma1=eye(d);
Sigma2=eye(d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting the sampled cold chain after MC3

t=length(Out);
plot(1:length(Out),Out);
ylabel('simulated value');
xlabel('iterate');

figure;
d=2;
m1=repmat(0,1,d);
m2=repmat(5,1,d;
Sigma1=eye(d);
Sigma2=eye(d);
p=0.5;

applyToGivenRow= @(func, matrix) @(row) func(matrix(row,:));
applyToRows=@ (func, matrix) arrayfun(applyToGivenRow(func,matrix),1:size(matrix,1))';
[x1 x2]= meshgrid(linspace(-3,7,100)',linspace(-3,7,100)');
X=[x1(:) x2(:)];
h=@(x) pdfmix(x,m1,m2,Sigma1,Sigma2,p);
Q=applyToRows(h,X);

Mixplot=surf(x1,x2,reshape(Q,100,100));
%imrotate(Mixplot,20);
xdir=[0 0 1];
rotate(Mixplot,xdir,15);
