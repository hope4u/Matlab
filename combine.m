clear;
for x=1:500
    load(sprintf('MIMOProj.%i.mat',x-1))
    R_ac_new(:,:,:,(x-1)*cdf+1:x*cdf) = R_ac;
    R_sum_new(:,:,:,(x-1)*cdf+1:x*cdf) = R_sum;
    R_max_new(:,:,:,(x-1)*cdf+1:x*cdf) = R_max;
    R_min_new(:,:,:,(x-1)*cdf+1:x*cdf) = R_min;
    P_tot_new(:,:,:,(x-1)*cdf+1:x*cdf) = P_tot;
    Error_new(:,:,(x-1)*cdf+1:x*cdf) = Error;
    H_err_new(:,:,:,(x-1)*cdf+1:x*cdf) = H_err;
end
R_ac = R_ac_new;
R_sum = R_sum_new;
R_max = R_max_new;
R_min = R_min_new;
P_tot = P_tot_new;
Error = Error_new;
H_err = H_err_new;

clear R_ac_new;
clear R_sum_new;
clear R_max_new;
clear R_min_new;
clear P_tot_new;
clear Error_new;
clear H_err_new;
clear x;

save combined