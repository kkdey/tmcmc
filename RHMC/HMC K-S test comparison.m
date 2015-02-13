L=100; % No. of runs for K-S comparison
HMC=cell(2,L);

% Assignment of the values
addpath('C:\Users\Kushal\Documents\MATLAB\RHMC');

d=50; % Dimensionality of the problem
nsamples=5000; % No. of iterations in each step
mu=zeros(1,d);
x0=repmat(5,1,d); % the starting point of position 
Sigma=eye(d);
M=30; % The number of leapfrog steps used 


% HMC run

for l=1:L
    smpl=zeros(nsamples,d);
    smpl(1,:)=x0;
    i=2;
    counter1=0;
    while i<=nsamples
        smpl(i,:)=HMC1(M,smpl((i-1),:),mu,Sigma);
        if smpl(i,:)==smpl((i-1),:)
            counter1=counter1+1;
        end
        i=i+1;
    end
    HMC{1,l}=smpl;
end

% RHMC run

for l=1:L
    smpl_RHMC=zeros(nsamples,d);
    smpl(1,:)=x0;
    i=2;
    counter2=0;
    while i<=nsamples
        smpl_RHMC(i,:)=RHMC1(M,smpl((i-1),:),mu,Sigma);
        if smpl_RHMC(i,:)==smpl_RHMC((i-1),:)
            counter2=counter2+1;
        end
        i=i+1;
    end
    HMC{2,l}=smpl_RHMC;
end


plot(1:nsamples,HMC{1,1}(:,2),'b');
hold on;
plot(1:nsamples,HMC{2,1}(:,2),'r');
hold off; 

KS_HMC=zeros(size(HMC{1,1},2),size(HMC{1,1},1));
KS_RHMC=zeros(size(HMC{1,1},2),size(HMC{1,1},1));


for n=1:size(HMC{1,1},2)
    for m=1:size(HMC{1,1},1)
        
        w=zeros(1,L);
        v=zeros(1,L);
        for l=1:L
            w(l) = HMC{1,l}(m,n);
            v(l) = HMC{2,l}(m,n);
        end
        
        
        samp = normrnd(mu(n),Sigma(n,n),1,L);
        [~,~,ks1]= kstest2(w,samp);
        [~,~,ks2]= kstest2(v,samp);
        KS_HMC(n,m)=ks1;
        KS_RHMC(n,m)=ks2;
    end
end

plot(1:nsamples,KS_HMC(1,:),'b');
hold on;
plot(1:nsamples,KS_RHMC(1,:),'r');
hold off;   
    
prop=zeros(1,d);
prop1=zeros(1,d);
nsamples1=nsamples/4;
for j=1:d
    prop(j)=sum(KS_HMC(j,:)>KS_RHMC(j,:))/nsamples;
    prop1(j)=sum(KS_HMC(j,nsamples1:nsamples)>KS_RHMC(j,nsamples1:nsamples))/(nsamples-nsamples1);
end

sum1=sum(prop>0.5);
sum2=sum(prop1>0.5);
