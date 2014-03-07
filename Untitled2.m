%% channel settings
N = 20;M = 20;
P = eye(N); % PowerMatrix

SNR = -20:5:20;
SNRLinear = 10^(SNR./10);
sigma = trace(P)./(N.*SNR);

Type={'LMMSE';'MMSE_VBLAST'};
Type=Type(2);



        %% Calculations
        SINR
        Rate = real(log2(SINR+1))
        R = sum(Rate)
%         R_ac = real(log2(det(eye(N)+sqrtm(P)*H'*sigma(j)^(-1)*H*sqrtm(P))))
        R_ac = real(log2(det(Phi)))
        
        
        

        %% Ploting
        
        figure(1)
        clf
        hold on
        subplot(1,2,1);
        hold on
        plot(SNR(j),R_ac_wf,'b*');
        plot(SNR(j),R_ac,'r.');
        subplot(1,2,2);
        hold on
        plot(SNR(j),R_ac_wf-R_ac,'go');


