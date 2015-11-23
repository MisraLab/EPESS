data {
int number_mixtures;
int dimension;
int dimension_sq;
real mixture_weights[number_mixtures];
matrix[number_mixtures, dimension] mixture_means;
real mix_cov[dimension_sq, number_mixtures];
}



parameters {
vector[dimension] x; 
}

model {
   real ps[number_mixtures];
   matrix[dimension, dimension] sigma;

   for (k in 1:number_mixtures) {
        for (i in 1:dimension){
             for (j in 1:dimension){
                 sigma[i,j] <- mix_cov[i+ (j-1)*dimension, k];
            }
       }
       
       ps[k] <- log(mixture_weights[k])
                 + multi_normal_log(x, mixture_means[k] , sigma);
   }
   increment_log_prob(log_sum_exp(ps));
   }