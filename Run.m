%% channel settings
clear;
close all;
seed_start=10;
randn('state',seed_start);
N = 4;M = 4;
P = eye(N); % PowerMatrix

SNR = linspace(15,55,40);
% SNR = 0;
SNRLinear = 10.^(SNR./10);

Type={'LMMSE'};
%Type:      receiver type
%     'LMMSE'               Linear MMSE equalizer
%     'MMSE_VBLAST'         MMSE with SIC (optimal receiver)
Optimizer={'fodorPrecoding'};
%Optimizer: 
%     'none'                no Power optimization
%     'wf'                  waterfilling and SVD precoding
%     'grad'                gradient Search for sumPower constraint
%     'sp_iwf'              sumPower constraint waterfilling
%     'sp_iwf_paper'        jindal's sumPower waterfilling
%     'fodor'               fodor's aproach with fairness constraints
%     'fodorPrecoding'      fodor's Precoding Optimization

%% run
[SINR, Phi] = MIMO_Transceiver(M,N,P,SNR,Type,Optimizer);

%% calculate Rate
R_ch = zeros(N,length(SNR),length(Type),length(Optimizer));
R_sum = zeros(length(SNR),length(Type),length(Optimizer));
R_ac = zeros(length(SNR),length(Type),length(Optimizer));

for j=1:length(SNR) %iterate over SNR
    for i = 1:length(Type) %iterate over Type   
        for k = 1:length(Optimizer) %iterate over Optimizer
            R_ch(:,j,i,k) = real(log2(SINR(:,j,i,k)+1));
            R_sum(j,i,k) = sum(R_ch(:,j,i,k));
            R_ac(j,i,k) = real(log2(det(Phi(:,:,j,i,k))));
        end
    end
end

%% plot

% Achievable Rate
% if exist('fRate','var'); figure(fRate); else fRate = figure; end; clf;
figure
hold on
for k=1:length(Optimizer)
    plot(SNR(:),R_sum(:,1,k),'color',[(length(Optimizer)-k)/(length(Optimizer)) k/(length(Optimizer)) 0]);
end
legend(Optimizer);
hold off