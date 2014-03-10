%% channel settings
N = 10;M = 20;
P = eye(N); % PowerMatrix

SNR = linspace(-20,0,10);
SNRLinear = 10.^(SNR./10);

Type={'LMMSE';'MMSE_VBLAST'};
Optimizer={'none';'wf';'ra_wf'};

% Optimizer=Optimizer{1};
% Optimizer=Optimizer{2};
Optimizer=Optimizer([1 3]);


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
figure(1);clf
hold on
plot(SNR(:),R_ac(:,1,1),'r');
plot(SNR(:),R_ac(:,1,2),'b');

hold off