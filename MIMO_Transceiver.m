function [SINR,Phi,P_ret] = MIMO_Transceiver(M,N,P,SNR,Type,Optimizer)

SNRLinear = 10.^(SNR./10);
sigma = trace(P)./(N.*SNRLinear);

%% Channel
[H] = MIMO_Channel(M,N,sigma(1));


%% Decoder
SINR = zeros(N,length(SNR),length(Type),length(Optimizer));
Phi = zeros(N,N,length(SNR),length(Type),length(Optimizer));
P_ret = zeros(N,N,length(SNR),length(Type),length(Optimizer));


for j=1:length(SNR) %iterate over SNR   
    for k = 1:length(Optimizer) %iterate over Optimizer
        switch Optimizer{k}
            case 'none'
                H_op = H;
                P_op = P;
            case 'wf'
                [ H_op,P_op ] = waterFilling(H,P,sigma(j));
            case 'sp_iwf'
                [ H_op,P_op ] = sumPower_iterativeWaterFilling( H,P,sigma(j) );
            case 'sp_iwf_paper'
                H_in = zeros(M,1,N);
                for n=1:N
                    H_in(:,:,n) = H(:,n);
                end
                [ ~,P_out ] = iterative_waterfill(H_in,trace(P),50);
                P_op = P;
                for n=1:N
                    P_op(n,n)=P_out(:,:,n);
                end
                H_op = H;
            case 'numericalGrad'
                [P_op, gradient] = numericalGradient(H,P,sigma(j));
                H_op = H;
                
            case 'analyticalGrad'
                [P_op, gradient] = analyticalGradient(H,P,sigma(j));
                H_op = H;
                
            case 'fodor'
                P_op = fodorloops(H,P,sigma(j));
                H_op = H;
                
            case 'fodorPrecoding'
                P_op = fodorPrecodingOptimization(H,P,sigma(j));
                H_op = H;
                
            case 'fodorPrecoding2'
                P_op = fodorPrecodingOptimization2(H,P,sigma(j));
                H_op = H;
                
            case 'minmax'
                [P_op, gradient, diffToTgt] = minmaxSINR(H,P,sigma(j));
%                 [P_op2, gradient2, diffToTgt2] = minPowerAnalytical( H,P,sigma(j) );
%                 P_op3 = minPowerClosedForm( H,P,sigma );
                H_op = H;
                
%             case 'minPower'
%                 P_op = minPower_rateConst(H,P,sigma(j));
%                 H_op = H;
        end
    
        for i = 1:length(Type) %iterate over Type               
            [Phi(:,:,j,i,k),SINR(:,j,i,k)] = MIMO_Receiver(H_op,P_op,sigma(j),Type{i});
            P_ret(:,:,j,i,k) = P_op;
        end
    end
end
