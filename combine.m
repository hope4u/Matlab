clear;
cdf = 10;
for i=1:100
    load(sprintf('MIMOProj.%i.mat',i-1))
    R_ac_new(:,:,:,(i-1)*cdf+1:i*cdf) = R_ac;
    R_ch_new(:,:,:,:,(i-1)*cdf+1:i*cdf) = R_ch;
    R_sum_new(:,:,:,(i-1)*cdf+1:i*cdf) = R_sum;
    H_all_new(:,:,:,:,:,(i-1)*cdf+1:i*cdf) = H_all;
    P_all_new(:,:,:,:,:,(i-1)*cdf+1:i*cdf) = P_all;
    Phi_all_new(:,:,:,:,:,(i-1)*cdf+1:i*cdf) = Phi_all;
    SINR_all_new(:,:,:,:,(i-1)*cdf+1:i*cdf) = SINR_all;
end
R_ac = R_ac_new;
R_ch = R_ch_new;
R_sum = R_sum_new;
H_all = H_all_new;
P_all = P_all_new;
Phi_all = Phi_all_new;
SINR_all = SINR_all_new;

clear R_ac_new;
clear R_ch_new;
clear R_sum_new;
clear H_all_new;
clear P_all_new;
clear Phi_all_new;
clear SINR_all_new;