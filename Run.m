%% channel settings
N = 10;M = 30;
P = eye(N); % PowerMatrix

SNR = linspace(-15,5,40);
SNRLinear = 10.^(SNR./10);

Type={'LMMSE';'MMSE_VBLAST'};
Optimizer={'none';'sp_iwf';'grad'};  %Optimizer={'none';'wf';'grad';'sp_iwf';'sp_iwf_paper'};


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
if exist('fRate','var'); figure(fRate); else fRate = figure; end; clf;
hold on
plot(SNR(:),R_ac(:,1,1),'r');
plot(SNR(:),R_ac(:,1,2),'g');
plot(SNR(:),R_ac(:,1,3),'b');
hold off