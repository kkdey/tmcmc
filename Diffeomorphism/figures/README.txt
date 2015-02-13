
This is a README file for the figures I have sent you..


The figures labeled "KS_compare_cauchy_dist_d_50_batch.. " correspond to the KS test for the first 20,000 iterations

for one particular dimension of a 50 dimensional distribution whose individual components are independent cauchy with 

mean 0 and scale 1. (batch number varies from 1 to 3 or 1 to 4 depending on how many replicates I have run, they are 

basically the same thing....just choose one among these figures)



The figures labeled "KS_compare_t_dist_d_50_batch.. " correspond to the KS test for the first 20,000 iterations

for one particular dimension of a 50 dimensional distribution whose individual components are independent t with degrees

of freedom equal to 10 and centered at 0. (batch number varies from 1 to 3 or 1 to 4 depending on how many replicates I

have run, they are basically the same thing....just choose one among these figures)


The figures labeled "KS_compare_t_dist_d_10_batch.. " correspond to first few iterations for the same thing as above but 

for 10 such independent t components, this is a 10 dimensional distribution.


The figures labeled "KS_MCMC_cauchy_w_and_wo_diffeomorphism_d_50_batch.." is the KS plots for the RWMH runs with and 

and without diffeomorphism for the Cauchy(0,1) distribution for dimension 50.Similarly 

"KS_MCMC_t_dist_w_and_wo_diffeomorphism_d_50_batch.." is the comparison with and without diffeomorphism for the t 

distribution centered at 0 and degrees of freedom 10.

The figures labeled "KS_TMCMC_cauchy_w_and_wo_diffeomorphism_d_50_batch.." is the KS plots for the TMCMC runs with and 

and without diffeomorphism for the Cauchy(0,1) distribution for dimension 50.Similarly 

"KS_TMCMC_t_dist_w_and_wo_diffeomorphism_d_50_batch.." is the comparison with and without diffeomorphism for the t 

distribution centered at 0 and degrees of freedom 10.


The figures labeled "KS_Mult_.._RWMH_TMCMC_compare_d_50_batch..." has two alternative distributions (t and cauchy)..in 

both cases, the mean/location is assumed to be 0, the degrees of freedom for t distribution to be 10, and the variance 

covariance matrix is assumed to be of the form

Sigma= 0.7 I + 0.3 1 1^{T}  

What we compare in these images are the RWMH and the TMCMC chains under diffeomorphism for these distributions.

"KS_TMCMC_Mult_..._w_and_wo_diffeomorphism_d_50_batch_..." represent the comparison of the chains with and without 

diffeomorphism for the Cauchy and the t distribution as mentioned above. 

"KS_MCMC_Mult_..._w_and_wo_diffeomorphism_d_50_batch_..." is the analogous representation for the RWMH chains.




