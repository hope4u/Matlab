function [SINR,Phi] = MIMO_Transceiver(M,N,P,SNR,Type,Optimizer)

SNRLinear = 10.^(SNR./10);
sigma = trace(P)./(N.*SNRLinear);

%% Channel
[H] = MIMO_Channel(M,N,sigma(1));


%% Decoder
SINR = zeros(N,length(SNR),length(Type),length(Optimizer));
Phi = zeros(N,N,length(SNR),length(Type),length(Optimizer));

for k = 1:length(Optimizer) %iterate over Optimizer
    for j=1:length(SNR) %iterate over SNR   
        switch Optimizer{k}
            case 'none'
                H_op = H;
                P_op = P;
            case 'wf'
                [ H_op,P_op ] = waterFilling(H,P,sigma(j));
            case 'grad'
                P_op = gradient(H,P,sigma(j));
                H_op = H;
        end
    
        for i = 1:length(Type) %iterate over Type               
            [Phi(:,:,j,i,k),SINR(:,j,i,k)] = MIMO_Receiver(H_op,P_op,sigma(j),Type{i});
        end
    end
end
