Haario=cell(2,L);

% Assignment of the values

d=2;
L=100;
mu=zeros(1,d);
mu1=repmat(0,1,d);
s_d=((2.42)^(2))/d;
x0=repmat(10,1,d);
epsilon1=0.10;
Sigma=eye(d);
Sigma1=eye(d);
nsamples=1000;
% The function
pdf=@(x)mvnpdf(x,mu,Sigma); % pdf of target density
proppdf=@(x,y)mvnpdf(x,mu1,Sigma1); % proposal distribution pdf 
proprnd = @(x)mvnrnd(mu1,Sigma1); % random number generator prop density

% Haario MCMC run

for l=1:L
    smpl=zeros(nsamples,d);
    smpl(1,:)=x0;
    i=2;
    counter=0;
    while i<=nsamples
        epsilon= mvnrnd(mu1, Sigma1,1);
        y= smpl((i-1),:)+epsilon;
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
        
        Sigma1= cov(smpl(1:(i-1),:))*(s_d)+ epsilon1*s_d*eye(d);
        proppdf=@(x,y)mvnpdf(x,mu1,Sigma1); % proposal distribution pdf
        proprnd = @(x)mvnrnd(mu1,Sigma1);
    end
    
    Haario{1,l}=smpl;
end


d=2;
L=100;
mu=zeros(1,d);
mu1=zeros(1,d);
s_d=((2.42)^(2))/d;
x0=repmat(10,1,d);
epsilon1=0.10;
Sigma=eye(d);
Sigma1=eye(d);
nsamples=1000;
% The function
pdf=@(x)mvnpdf(x,mu,Sigma); % pdf of target density
proppdf=@(x,y)mvnpdf(x,mu1,Sigma1); % proposal distribution pdf 
proprnd = @(x)mvnrnd(mu1,Sigma1);


% Haario TMCMC run



for l=1:L
    smpl_TMCMC=zeros(nsamples,d);
    x0=repmat(10,1,d);
    smpl_TMCMC(1,:)=x0;
    s_d=((2.42)^(2))/d;
    i=2;
    counter=0;
    store=zeros(1,nsamples);
    mu1=zeros(1,d);
    Sigma=eye(d);
    pdf=@(x)mvnpdf(x,mu,Sigma);
    Sigma1=eye(d);
    
    sigmaprop=Sigma1(1,1);
    epsilon1=0.10;

    while i<=nsamples
        epsilon= abs(normrnd(0,sigmaprop));
        store(i)=epsilon;
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
            sigmaprop=(s_d)*var(store(2:i))+epsilon1*(s_d);
            i=i+1;
            
        else
            smpl_TMCMC(i,:)=smpl_TMCMC((i-1),:);
            sigmaprop=(s_d)*var(store(2:i))+epsilon1*(s_d);
            i=i+1;
            
        end
    end
    
    Haario{2,l}=smpl_TMCMC;
end

plot(1:nsamples,Haario{1,1}(:,2),'b');
hold on;
plot(1:nsamples,Haario{2,1}(:,2),'r');
hold off; 

KS_MCMC=zeros(size(Haario{1,1},2),size(Haario{1,1},1));
KS_TMCMC=zeros(size(Haario{1,1},2),size(Haario{1,1},1));


for n=1:size(Haario{1,1},2)
    for m=1:size(Haario{1,1},1)
        
        w=zeros(1,L);
        v=zeros(1,L);
        for l=1:L
            w(l) = Haario{1,l}(m,n);
            v(l) = Haario{2,l}(m,n);
        end
        
        
        samp = normrnd(mu(n),Sigma(n,n),1,L);
        [~,~,ks1]= kstest2(w,samp);
        [~,~,ks2]= kstest2(v,samp);
        KS_MCMC(n,m)=ks1;
        KS_TMCMC(n,m)=ks2;
    end
end

plot(1:nsamples,KS_MCMC(1,:),'b');
hold on;
plot(1:nsamples,KS_TMCMC(1,:),'r');
hold off;   
    
prop=zeros(1,d);
prop1=zeros(1,d);
nsamples1=nsamples/4;
for j=1:d
    prop(j)=sum(KS_MCMC(j,:)>KS_TMCMC(j,:))/nsamples;
    prop1(j)=sum(KS_MCMC(j,nsamples1:nsamples)>KS_TMCMC(j,nsamples1:nsamples))/(nsamples-nsamples1);
end

sum1=sum(prop>0.5);
sum2=sum(prop1>0.5);