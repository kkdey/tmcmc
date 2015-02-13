% Random Dive Metropolis Hastings Algorithm (TMCMC approach)- with additive
% TMCMC as within a neighbourhood of 0

% Assignment of the variables

d=30;
L=30;
propvar=power(6,2)/d;
x0=repmat(1,1,d);


% Algorithm construction

Dive=cell(2,L);
mu=zeros(1,d);
mu1=repmat(0.3,1,d);
mu2=repmat(0,1,d);


Sigma=eye(d);
Sigma1=propvar*eye(d);
nsamples=40000;
% The function
pdf=@(x)mvnpdf(x,mu,Sigma); % pdf of target density
proppdf=@(x,y)mvnpdf(x,mu1,Sigma1); % proposal distribution pdf 
proprnd = @(x)mvnrnd(mu1,Sigma1); % random number generator prop density

%MCMC run
for l=1:L
    smpl=zeros(nsamples,d);
    smpl(1,:)=x0;
    i=2;
    counter=0;
    while i<=nsamples
        epsilon= mvnrnd(mu2, Sigma1,1);
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
            continue;
        end
    end
    
    Dive{1,l}=smpl;
end
 acceptrateMCMC=counter/nsamples;
 fprintf(' %d,  Acceptance ratio of process \n', acceptrateMCMC);
 
 
% Random Dive Metropolis Hastings algorithm
propvar2=power(0.3,2)/d;
Sigma1=propvar2*eye(d);
sigmaprop=Sigma1(1,1);
sigmaprop2=sqrt(power(6,2)/d);

for l=1:L
    smpl_TMCMC=zeros(nsamples,d);
    smpl_TMCMC(1,:)=x0;
    i=2;
    counter=0;
    
    
    while i<=nsamples
        
        v=rand(1);
        
        if v<0.5
            
            epsilon=2;
            while epsilon >0.90 ||epsilon < 0.05
                epsilon= abs(normrnd(mu1(1),sqrt(sigmaprop)));
            end
            
            
            
            y=zeros(1,d);
            counter_array=zeros(1,4);
            ep=[epsilon -epsilon 1/epsilon -1/epsilon ];
            for x=1:d
                u=rand(1);
                if u<0.25
                    y(x)=smpl_TMCMC((i-1),x)*ep(1);
                    counter_array(1)= counter_array(1)+1;
                else if 0.25 <= u < 0.5
                        y(x)=smpl_TMCMC((i-1),x)*ep(2);
                        counter_array(2)= counter_array(2)+1;
                    else if 0.5 <= u < 0.75
                            y(x)=smpl_TMCMC((i-1),x)*ep(3);
                            counter_array(3)= counter_array(3)+1;
                        else
                            y(x)=smpl_TMCMC((i-1),x)*ep(4);
                            counter_array(4)= counter_array(4)+1;
                        end
                    end
                end
            end
            
            accept_rate= min(1,(pdf(y)*power(ep(1),counter_array(1)+counter_array(2))* power(ep(3),counter_array(3)+counter_array(4))/pdf(smpl_TMCMC((i-1),:))));
            
        end
        
        if v>0.5
            y=zeros(1,d);
            epsilon= abs(normrnd(0,sigmaprop2));
            for x=1:d
                u=rand(1);
                if u<0.5
                    y(x)=smpl_TMCMC((i-1),x)+epsilon;
                else
                    y(x)=smpl_TMCMC((i-1),x)-epsilon;
                end
            end
            
            accept_rate= min(1,(pdf(y)/pdf(smpl_TMCMC((i-1),:))));
        end
        
        
        
        
        
        w=rand(1);
        if w<accept_rate
            smpl_TMCMC(i,:)=y;
            i=i+1;
            counter=counter+1;
        else
            smpl_TMCMC(i,:)=smpl_TMCMC((i-1),:);
            i=i+1;
            continue;
        end
    end
    
    Dive{2,l}=smpl_TMCMC;
end

plot(1:nsamples,Dive{1,1}(1:nsamples,1),'b');
hold on;
plot(1:nsamples,Dive{2,1}(1:nsamples,1),'r');
hold off; 
legend('RWMH','RDMH');
xlabel('Iterate');
ylabel('value');



KS_MCMC=zeros(size(Dive{1,1},2),size(Dive{1,1},1));
KS_TMCMC=zeros(size(Dive{1,1},2),size(Dive{1,1},1));


for n=1:size(Dive{2,1},2)
    for m=1:size(Dive{2,1},1)
        
        w=zeros(1,L);
        v=zeros(1,L);
        for l=1:L
            w(l) = Dive{1,l}(m,n);
            v(l) = Dive{2,l}(m,n);
        end
        
        
        samp = normrnd(mu(n),Sigma(n,n),1,L);
        [~,~,ks1]= kstest2(w,samp);
        [~,~,ks2]= kstest2(v,samp);
        KS_MCMC(n,m)=ks1-0.1;
        KS_TMCMC(n,m)=ks2-0.1;
    end
end

plot(1:nsamples,KS_MCMC(10,:),'b');
hold on;
plot(1:nsamples,KS_TMCMC(5,:),'r');
hold off;   
    
m=4;
n=1000;



plot(1:n,KS_MCMC(m,1:n)+repmat(0.02,1,n),'b');
hold on;
plot(1:n,KS_TMCMC(m,1:n)+repmat(0.07,1,n),'r');
hold off;  
legend('RWMH','ATMCMC');
xlabel('Iterate');
ylabel('K-S statistic');
ylim([0 1]);

prop=zeros(1,d);
prop1=zeros(1,d);
nsamples1=nsamples/4;
for j=1:d
    prop(j)=sum(KS_MCMC(j,:)>KS_TMCMC(j,:))/nsamples;
    prop1(j)=sum(KS_MCMC(j,nsamples1:nsamples)>KS_TMCMC(j,nsamples1:nsamples))/(nsamples-nsamples1);
end

sum1=sum(prop>0.5);
sum2=sum(prop1>0.5);

KS_M_avg= mean(mean(KS_MCMC));
KS_T_avg= mean(mean(KS_TMCMC));

acceptrateTMCMC=counter/nsamples;
fprintf('%d,  Acceptance ratio of process \n', acceptrateTMCMC);
fprintf('%d, Avg K-S test value for MCMC \n', KS_M_avg);
fprintf('%d, Avg K-S test value for RDMH \n',KS_T_avg);
