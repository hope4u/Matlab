%% channel settings
clear;
close all;
seed_start=10;
randn('state',seed_start);
N = 12;M = 12;
P = eye(N); % PowerMatrix

SNR = 20;   
SNRLinear = 10.^(SNR./10);
cdf=10000;


Type={'LMMSE'};
%Type:      receiver type
%     'LMMSE'               Linear MMSE equalizer
%     'MMSE_VBLAST'         MMSE with SIC (optimal receiver)
Optimizer={'none'};
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

for j=1:cdf
%% run
[SINR, Phi] = MIMO_Transceiver(M,N,P,SNR,Type,Optimizer);

%% calculate Rate


    for i = 1:length(Type) %iterate over Type   
        for k = 1:length(Optimizer) %iterate over Optimizer
            R_ch(:,j,i,k) = real(log2(SINR(:,1,i,k)+1));
            R_sum(j,i,k) = sum(R_ch(:,j,i,k));
            R_ac(j,i,k) = real(log2(det(Phi(:,:,1,i,k))));
        end
    end
end

%% plot

% Achievable Rate
ecdf(min(R_ch))
figure
hold on
for k=1:length(Optimizer)
    plot(SNR(:),R_sum(:,1,k),'color',[1-(k-1)/(length(Optimizer)-1) (k-1)/(length(Optimizer)-1) 0]);
end
hold off
