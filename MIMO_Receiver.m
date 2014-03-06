function [ x_rec,Phi,SINR ] = MIMO_Receiver( y,H,P,sigma, Type )
M = size(H,1);
N = size(H,2);

switch Type
    case 'LMMSE'
        if size(sigma) == [M,M]
            K_n = sigma;
        elseif size(sigma,1) == M
            K_n = diag(sigma);
        else
            K_n = eye(M)*sigma;
        end
        
        Phi = (P*H'*K_n^(-1)*H+eye(N));
        G = Phi^(-1) *H'*K_n^(-1);
        SINR = 1./diag(Phi^(-1)) - 1;
        x_rec = G*y;
        
    case 'MMSE_VBLAST'
        if size(sigma) == [M,M]
            K_n = sigma;
        elseif size(sigma,1) == M
            K_n = diag(sigma);
        else
            K_n = eye(M)*sigma;
        end
        
        Phi = (P*H'*K_n^(-1)*H+eye(N));
        G = Phi^(-1) *H'*K_n^(-1);
        x_rec = G*y;

        SINR = zeros(N,1);
        index = [1:N]';
        for i = 1:N
            currentPhi = (P*H'*K_n^(-1)*H+eye(size(index,1)));

            currentSINR = 1./diag(currentPhi^(-1))-1;
            [maxSINR, maxSINR_index] = max(currentSINR);
            
            SINR_index = index(maxSINR_index);
            SINR(SINR_index) = maxSINR;
           
            index(maxSINR_index)=[];
            H(:,maxSINR_index)=[];
            P(:,maxSINR_index)=[]; P(maxSINR_index,:)=[];
        end
        
end
end

