%% channel settings
N = 20;M = 20;
p = 1+0*randn(N,1);
P = diag(p); % PowerMatrix

SNR = .1:.1:3;
sigma = trace(P)./(N.*SNR);

Type={'LMMSE';'MMSE_VBLAST'};
Type=Type(2);


numOfSymb = N;
bitsPerSymb = 2;
numOfBits = 2*N;

%% generate Data
data = randi([0,1],numOfBits,1)*2-1;

x = reshape(data,[],2) * [1;1i] * 1/sqrt(2);

%% send Data
[y,H] = MIMO_Channel(x,M,N,sigma(1));

figure(1)
clf
hold on
for j=1:length(SNR)
    %% receive Data
    for i = 1:size(Type)
        fprintf(1,'\n%s:\n',Type{i});
        fprintf(1,'%s\n',repmat('-',1,length(Type{i})+1));
        [x_rec,Phi,SINR] = MIMO_Receiver(y,H,P,sigma(j),Type{i});
        e = x_rec-x;

        %% Calculations
        SINR
        Rate = real(log2(SINR+1))
        R = sum(Rate)
        R_ac = real(log2(det(eye(N)+sqrtm(P)*H'*sigma(j)^(-1)*H*sqrtm(P))))
%         R_ac = real(log2(det(Phi)))

        %% Water-Filling
        [U,S,V]=svd(H);
        s = diag(S'*S);
        v = trace(P);
        for k=1:N
            if v > sigma(j)/s(k)
                current_v = (sum(sigma(j)./s(1:k)) + trace(P))/k;
                if current_v > sigma(j)/s(k)
                    v = current_v;
                end
            end
        end
        P_wf = v - sigma(j)./(S'*S);
        P_wf(P_wf<0)=0;
        P_wf=V*P_wf*V';
        H_wf=H;


        fprintf(1,'\nWith Water-Filling:\n');
        fprintf(1,'\n%s:\n',Type{i});
        fprintf(1,'%s\n',repmat('-',1,length(Type{i})+1));
        [x_rec_wf,Phi_wf,SINR_wf] = MIMO_Receiver(y,H,P_wf,sigma(j),Type{i});
        e = x_rec_wf-x;

        %% Calculations
        SINR_wf
        Rate_wf = real(log2(SINR_wf+1))
        R_wf = sum(Rate_wf)
        R_ac_wf = real(log2(det(eye(N)+sqrtm(P_wf)*H_wf'*sigma(j)^(-1)*H_wf*sqrtm(P_wf))))
%         R_ac_wf = real(log2(det(Phi)))
        subplot(1,2,1);
        hold on
        plot(SNR(j),R_ac_wf,'b*');
        plot(SNR(j),R_ac,'r.');
        subplot(1,2,2);
        hold on
        plot(SNR(j),R_ac_wf-R_ac,'go');
    end
end
hold off
%% Postprocessing
