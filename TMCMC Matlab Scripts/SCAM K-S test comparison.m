% SCAM K-S test comparison
L=100;
SCAM=cell(2,L);


% Assignment of the values

d=3;
L=100;
mu=zeros(1,d);
mu1=zeros(1,d);
x0=repmat(10,1,d);
Sigma=eye(d);
Sigma1=eye(d);
nsamples=1000;
% The function
pdf=@(x)mvnpdf(x,mu,Sigma); % pdf of target density
proppdf=@(x,y)mvnpdf(x,mu1,Sigma1); % proposal distribution pdf 
proprnd = @(x)mvnrnd(mu1,Sigma1); % random number generator prop density

%SCAM MCMC run
for l=1:L
    smpl=zeros(nsamples,d);
    smpl(1,:)=x0;
    
    i=2;
    counter=0;
    Variance=zeros(1,d);
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
        
        if i<=10
            for x=1:d
                Variance(x)=25;
            end
        else
            for x=1:d
                Variance(x)=(power(2.4,2))*(var(smpl(1:(i-1),x))+0.05);
            end
        end
        
        Sigma1= diag(Variance,0);
        proppdf=@(x,y)mvnpdf(x,mu1,Sigma1); % proposal distribution pdf
        proprnd = @(x)mvnrnd(mu1,Sigma1); % random number generator prop density
     end
    SCAM{1,l}=smpl;
end

% SCAM TMCMC run

sigmaprop=Sigma1(1,1);

for l=1:L
    smpl_TMCMC=zeros(nsamples,d);
    smpl_TMCMC(1,:)=x0;
    i=2;
    counter=0;
    while i<=nsamples
        epsilon= abs(normrnd(0,sigmaprop));
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
            i=i+1;
        else
            smpl_TMCMC(i,:)=smpl_TMCMC((i-1),:);
            i=i+1;
        end
        
        
        var1=0;
        if i<=10
            var1=25;
        else
            for x=1:d
                var1=var1+(1/d)*(power(2.4,2))*(var(smpl_TMCMC(1:(i-1),x))+0.05);
            end
            
        end
        
        sigmaprop=var1;
        %     proppdf=@(x,y)mvnpdf(x,[0 0],Sigma1); % proposal distribution pdf
        %     proprnd = @(x)mvnrnd([0 0],Sigma1); % random number generator prop density
    end
    SCAM{2,l}=smpl_TMCMC;
end

plot(1:nsamples,SCAM{1,78}(:,1),'b');
hold on;
plot(1:nsamples,SCAM{2,78}(:,1),'r');
hold off; 
legend('RWMH path','TMCMC path');

KS_MCMC=zeros(size(SCAM{1,1},2),size(SCAM{1,1},1));
KS_TMCMC=zeros(size(SCAM{1,1},2),size(SCAM{1,1},1));


for n=1:size(SCAM{1,1},2)
    for m=1:size(SCAM{1,1},1)
        
        w=zeros(1,L);
        v=zeros(1,L);
        for l=1:L
            w(l) = SCAM{1,l}(m,n);
            v(l) = SCAM{2,l}(m,n);
        end
        
        
        samp = normrnd(mu(n),Sigma(n,n),1,L);
        [~,~,ks1]= kstest2(w,samp);
        [~,~,ks2]= kstest2(v,samp);
        KS_MCMC(n,m)=ks1;
        KS_TMCMC(n,m)=ks2;
    end
end

plot(1:nsamples,KS_MCMC(3,:),'b');
hold on;
plot(1:nsamples,KS_TMCMC(3,:),'r');
hold off;   
    
prop=zeros(1,d);
prop1=zeros(1,d);

% nsamples1=nsamples/2;
nsamples1=700;
for j=1:d
    prop(j)=sum(KS_MCMC(j,:)>KS_TMCMC(j,:))/nsamples;
    prop1(j)=sum(KS_MCMC(j,nsamples1:nsamples)>KS_TMCMC(j,nsamples1:nsamples))/(nsamples-nsamples1);
end

sum1=sum(prop>0.5);
sum2=sum(prop1>0.5);

    

    
