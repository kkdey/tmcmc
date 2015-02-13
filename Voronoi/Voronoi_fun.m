function [Out1 Out2] = Voronoi_fun(N,x0, Sigma1, sigmaprop, pdf  )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

W_MCMC=MCMC_fun(N,2,x0,pdf,Sigma1); % MCMC run 
W_TMCMC=TMCMC_fun(N,2,x0,pdf,sigmaprop); % TMCMC run

% figure; Voronoi_MCMC=voronoi(W_MCMC(:,1),W_MCMC(:,2));
% figure; Voronoi_TMCMC=voronoi(W_TMCMC(:,1),W_TMCMC(:,2));

[V1 C1]=voronoin(unique(W_MCMC,'rows'));
tess_area_MCMC=zeros(size(C1,1),1);



for i = 1 : size(C1,1)
    ind = C1{i}';
    tess_area_MCMC(i,1) = polyarea( V1(ind,1) , V1(ind,2) );
end

[V2 C2]=voronoin(unique(W_TMCMC,'rows'));
tess_area_TMCMC=zeros(size(C1,1),1);



for i = 1 : size(C2,1)
    ind = C2{i}';
    tess_area_TMCMC(i,1) = polyarea( V2(ind,1) , V2(ind,2) );
end

u_MCMC=unique(W_MCMC,'rows');
u_TMCMC=unique(W_TMCMC,'rows');

lab1=find(~isnan(tess_area_MCMC(:,1)));
lab2=find(~isnan(tess_area_TMCMC(:,1)));

%Netdata_MCMC=[u_MCMC(lab1,:) tess_area_MCMC(lab1,1)];
%Netdata_TMCMC=[u_TMCMC(lab2,:) tess_area_TMCMC(lab2,1)];

cell1=num2cell(u_MCMC(lab1,:),2);
cell2=num2cell(u_TMCMC(lab2,:),2);
tess_area_MCMC_lab=tess_area_MCMC(lab1,1);
tess_area_TMCMC_lab=tess_area_TMCMC(lab2,1);


pdfcell1= cell2mat(cellfun(@(x) pdf(x), cell1, 'UniformOutput', false));
pdfcell2= cell2mat(cellfun(@(x) pdf(x), cell2, 'UniformOutput', false));

sum_MCMC=0;
sum_TMCMC=0;

for i=1:length(pdfcell1)
    sum_MCMC=sum_MCMC+pdfcell1(i)*tess_area_MCMC_lab(i);
end

for i=1:length(pdfcell2)
  sum_TMCMC=sum_TMCMC+pdfcell2(i)*tess_area_TMCMC_lab(i);  

end

Out1=sum_MCMC;
Out2=sum_TMCMC;


end

