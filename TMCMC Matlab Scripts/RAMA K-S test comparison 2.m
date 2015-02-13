%RAMA K-S test comparison (2 partition)

RAMA=cell(2,L);

d=50;
L=100;
mu=zeros(1,d);
mu1=repmat(0,1,d);
x0=repmat(10,1,d);
Sigma=eye(d);
Sigma1=eye(d);
nsamples=1000;
% The function
pdf=@(x)mvnpdf(x,mu,Sigma); % pdf of target density
proppdf=@(x,y)mvnpdf(x,mu1,Sigma1); % proposal distribution pdf 
proprnd = @(x)mvnrnd(mu1,Sigma1); % random number generator prop density


for l=1:L
    D=20;
    J=100;
    smpRAMA=zeros(D,d);
    pdf=@(x)mvnpdf(x,mu,Sigma); % pdf of target density
    Sigma1=5*Sigma;
    mu1=zeros(1,d);
    proppdf=@(x,y)mvnpdf(x,mu1,Sigma1); % proposal distribution pdf
    proprnd = @(x)mvnrnd(mu1,Sigma1);
    delta=zeros(1,(J+1));
    a=zeros(1,(J+1));
    b=zeros(1,(J+1));
    a(1)=rand(1);
    b(1)=rand(1);
    SampleRAMA=zeros(1,d);
    
    for i=1:J
        counter1=0;
        counter2=0;
        Bota=zeros(1,D);
        smpRAMA(1,:)=SampleRAMA(end,:);
        for j=2:D
            y=smpRAMA((j-1),:);
            if (norm(y)<2)
                sigmaprop=exp(2*a(i));
                Bota(j)=1;
            else
                sigmaprop=exp(2*b(i));
                Bota(j)=2;
            end
            epsilon= normrnd(0,sigmaprop,1,d);
            %store(i)=epsilon;
            %sigmaprop=(s_d)*var(store(1:i))+epsilon*(s_d);
            y=y+epsilon;
            accept_rate= min(1,(pdf(y)/pdf(smpRAMA((j-1),:))));
            w=rand(1);
            if w<accept_rate
                smpRAMA(j,:)=y;
                if (Bota(j)==1)
                    counter1=counter1+1;
                else
                    counter2=counter2+1;
                end
                
            else
                smpRAMA(j,:)=smpRAMA((j-1),:);
                continue;
            end
        end
        alphahat1= counter1/length(find(Bota==1));
        alphahat2= counter2/length(find(Bota==2));
        delta(i)=min(0.01,sqrt(1/i));
        SampleRAMA=[SampleRAMA; smpRAMA(2:end,:)];
        if (alphahat1<0.439)
            a(i+1)=a(i)-delta(i);
        else
            a(i+1)=a(i)+delta(i);
        end
        
        if (alphahat2<0.439)
            b(i+1)=b(i)-delta(i);
        else
            b(i+1)=b(i)+delta(i);
        end
        
    end
    
    RAMA{1,l}=SampleRAMA;
end



for l=1:L
    D=20;
    J=100;
    smpRAMA=zeros(D,d);
    pdf=@(x)mvnpdf(x,mu,Sigma); % pdf of target density
    Sigma1=5*Sigma;
    mu1=zeros(1,d);
    proppdf=@(x,y)mvnpdf(x,mu1,Sigma1); % proposal distribution pdf
    proprnd = @(x)mvnrnd(mu1,Sigma1);
    delta=zeros(1,(J+1));
    a=zeros(1,(J+1));
    b=zeros(1,(J+1));
    a(1)=rand(1);
    b(1)=rand(1);
    SampleRAMA=zeros(1,d);
    
    for i=2:J
        counter1=0;
        counter2=0;
        Bota=zeros(1,D);
        smpRAMA(1,:)=SampleRAMA(end,:);
        for j=2:D
            y=smpRAMA((j-1),:);
            if (norm(y)<2)
                sigmaprop=exp(2*a(i));
                Bota(j)=1;
            else
                sigmaprop=exp(2*b(i));
                Bota(j)=2;
            end
            epsilon= abs(normrnd(0,sigmaprop));
            y=zeros(1,d);
            for x=1:d
                u=rand(1);
                if u<0.5
                    y(x)=smpRAMA((j-1),x)+epsilon;
                else
                    y(x)=smpRAMA((j-1),x)-epsilon;
                end
            end
            
            accept_rate= min(1,(pdf(y)/pdf(smpRAMA((j-1),:))));
            w=rand(1);
            if w<accept_rate
                smpRAMA(j,:)=y;
                if (Bota(j)==1)
                    counter1=counter1+1;
                else
                    counter2=counter2+1;
                end
                
            else
                smpRAMA(j,:)=smpRAMA((j-1),:);
                continue;
            end
        end
        alphahat1= counter1/length(find(Bota==1));
        alphahat2= counter2/length(find(Bota==2));
        delta(i)=min(0.01,sqrt(1/i));
        SampleRAMA=[SampleRAMA; smpRAMA(2:end,:)];
        if (alphahat1<0.439)
            a(i+1)=a(i)-delta(i);
        else
            a(i+1)=a(i)+delta(i);
        end
        
        if (alphahat2<0.439)
            b(i+1)=b(i)-delta(i);
        else
            b(i+1)=b(i)+delta(i);
        end
        
    end
    RAMA{2,l}=SampleRAMA;
end

plot(1:nsamples,RAMA{1,1}(:,48),'b');
hold on;
plot(1:nsamples,RAMA{2,1}(:,48),'r');
hold off; 

KS_MCMC=zeros(size(RAMA{1,1},2),size(RAMA{1,1},1));
KS_TMCMC=zeros(size(RAMA{1,1},2),size(RAMA{1,1},1));


for n=1:size(RAMA{1,1},2)
    for m=1:size(RAMA{1,1},1)
        
        w=zeros(1,L);
        v=zeros(1,L);
        for l=1:L
            w(l) = RAMA{1,l}(m,n);
            v(l) = RAMA{2,l}(m,n);
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

    
