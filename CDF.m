%% channel settings
clear;
close all;
seed_start=10;
randn('state',seed_start);
N = 4;M = 4;
P = eye(N); % PowerMatrix

SNR = 1;   
SNRLinear = 10.^(SNR./10);
cdf=1;


Type={'MMSE_VBLAST'};
%Type:      receiver type
%     'LMMSE'               Linear MMSE equalizer
%     'MMSE_VBLAST'         MMSE with SIC (optimal receiver)
Optimizer={'sp_iwf_paper','numericalGrad_VBlast'};
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

for j=1:cdf
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

%% plot

% Achievable Rate

figure
hold all
for k=1:length(Optimizer)
    ecdf(R_sum(:,1,k));%,'color',[(length(Optimizer)-k)/length(Optimizer) k/length(Optimizer) 0]);
end
hleg = legend(Optimizer);
set(hleg,'Interpreter','none')
set(hleg,'Location','Best')
hold off

% name = sprintf('%i%s',cdf,strjoin(Optimizer,'\b'));
save(filename)
