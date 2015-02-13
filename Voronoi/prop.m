W=5:5:100;
x0=[0 0];
Data=zeros(length(W),2);
for k=1:length(W)
    Sigma1=W(k)*eye(2);
    sigmaprop=W(k);
    [out1 out2]= Voronoi_fun(10000,x0,Sigma1, sigmaprop,pdf);
    Data(k,:)=[out1 out2];
end

figure; plot(W, Data(:,1), 'b'); xlabel('Proposal variance'); ylabel('Value of integral'); hold on;
plot(W, Data(:,2), 'r'); t=0:0.0001:1;
h4=plot([0 W],repmat(11.106,length(W)+1),'green-'); legend('MCMC','TMCMC', 'True'); ylim([10,30]);



    
