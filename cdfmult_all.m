function [] = cdfmult(filename,number)
%% channel settings
close all;
randn('state',number*10);
N = 12;M = 12;
P = eye(N); % PowerMatrix

SNR = [-10,0,10,20];
TgtRate = [1,2,10];
SNRLinear = 10.^(SNR./10); % target Rates for minpower
cdf=200;


Type={'LMMSE','MMSE_VBLAST'};
%Type:      receiver type
%     'LMMSE'               Linear MMSE equalizer
%     'MMSE_VBLAST'         MMSE with SIC (optimal receiver)
Optimizer={'none','wf','sp_iwf','sp_iwf_paper','numericalGrad_VBlast','numericalGrad','minmax','minmax_VBlast','powermin','fodorPrecoding'};
% Optimizer={'none','wf','sp_iwf','sp_iwf_paper','numericalGrad_VBlast'};%'numericalGrad','minmax'};

%Optimizer: 
%     'none'                no Power optimization
%     'wf'                  waterfilling and SVD precoding
%     'grad'                gradient Search for sumPower constraint
%     'sp_iwf'              sumPower constraint waterfilling
%     'sp_iwf_paper'        jindal's sumPower waterfilling
%     'fodor'               fodor's aproach with fairness constraints

R_ch = zeros(N,length(SNR),length(Type),length(Optimizer),cdf);
R_sum = zeros(length(SNR),length(Type),length(Optimizer),cdf);
R_ac = zeros(length(SNR),length(Type),length(Optimizer),cdf);
R_min = zeros(length(SNR),length(Type),length(Optimizer),cdf);
R_max = zeros(length(SNR),length(Type),length(Optimizer),cdf);
P_tot = zeros(length(SNR),length(Type),length(Optimizer),cdf);
Error = zeros(3,length(SNR),cdf);
H_err = zeros(N,N,length(SNR),cdf);

numSNR = length(SNR);
numType = length(Type);
numOptimizer = length(Optimizer);

for n=1:cdf
%% run
[SINR, Phi, P_ret, H] = MIMO_Transceiver(M,N,P,SNR,Type,Optimizer);

%% calculate Rate
% 
% P_all(:,:,:,:,:,n) = P_ret;
% H_all(:,:,n) = H;
% Phi_all(:,:,:,:,:,n) = Phi;
% SINR_all(:,:,:,:,n) = SINR;

    for j=1:numSNR
        for i = 1:numType %iterate over Type   
            for k = 1:numOptimizer %iterate over Optimizer
                R_ch(:,j,i,k,n) = real(log2(SINR(:,j,i,k)+1));
                R_sum(j,i,k,n) = sum(R_ch(:,j,i,k,n));
                R_ac(j,i,k,n) = real(log2(det(Phi(:,:,j,i,k))));
                
                if sum(strcmp({'fodorPrecoding2','powermin'},Optimizer{k}))
                    P_tot(j,i,k,n) = trace(P_ret(:,:,j,i,k));
                else
                    R_max(j,i,k,n) = max(R_ch(:,j,i,k,n));
                    R_min(j,i,k,n) = min(R_ch(:,j,i,k,n));
                end
                        

            end
        end
    end
    if R_sum(1,1,6,n)<8
        Error(:,1,n) = [number,1,n];
        H_err(:,:,1,n) = H;
        P_tot(1,1,6,n) = trace(P_ret(:,:,1,1,6));
    elseif R_sum(2,1,6,n)<20
        Error(:,2,n) = [number,2,n];
        H_err(:,:,2,n) = H;
        P_tot(2,1,6,n) = trace(P_ret(:,:,2,1,6));
    elseif R_sum(3,1,6,n)<40
        Error(:,3,n) = [number,3,n];
        H_err(:,:,3,n) = H;
        P_tot(3,1,6,n) = trace(P_ret(:,:,3,1,6));
    elseif R_sum(4,1,6,n)<60
        Error(:,4,n) = [number,4,n];
        H_err(:,:,4,n) = H;
        P_tot(4,1,6,n) = trace(P_ret(:,:,4,1,6));
    end
    
end

save(filename,'R_sum','R_ac','R_min','R_max','P_tot','SNR','Type','Optimizer','cdf','Error','H_err')

end
