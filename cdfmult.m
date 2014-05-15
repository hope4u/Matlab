function [] = cdfmult(filename,number)
%% channel settings
close all;

randn('state',number*10);
N = 12;M = 12;
P = eye(N); % PowerMatrix

SNR = 1;    
SNRLinear = 10.^(SNR./10);
cdf=10;


Type={'LMMSE'};
%Type:      receiver type
%     'LMMSE'               Linear MMSE equalizer
%     'MMSE_VBLAST'         MMSE with SIC (optimal receiver)
Optimizer={'numericalGrad'};
%Optimizer: 
%     'none'                no Power optimization
%     'wf'                  waterfilling and SVD precoding
%     'grad'                gradient Search for sumPower constraint
%     'sp_iwf'              sumPower constraint waterfilling
%     'sp_iwf_paper'        jindal's sumPower waterfilling
%     'fodor'               fodor's aproach with fairness constraints

R_ch = zeros(N,cdf,length(Type),length(Optimizer));
R_sum = zeros(cdf,length(Type),length(Optimizer));
R_ac = zeros(cdf,length(Type),length(Optimizer));

numType = length(Type);
numOptimizer = length(Optimizer);

parfor j=1:cdf
%% run
[SINR, Phi] = MIMO_Transceiver(M,N,P,SNR,Type,Optimizer);

%% calculate Rate


    for i = 1:numType %iterate over Type   
        for k = 1:numOptimizer %iterate over Optimizer
            R_ch(:,j,i,k) = real(log2(SINR(:,1,i,k)+1));
            R_sum(j,i,k) = sum(R_ch(:,j,i,k));
            R_ac(j,i,k) = real(log2(det(Phi(:,:,1,i,k))));
        end
    end
end

save(filename)

end
