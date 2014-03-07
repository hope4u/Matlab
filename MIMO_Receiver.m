function [Phi,SINR ] = MIMO_Receiver( H,P,sigma, Type )
M = size(H,1);
N = size(H,2);

if size(sigma) == [M,M]
    K_n = sigma;
elseif size(sigma,1) == M
    K_n = diag(sigma);
else
    K_n = eye(M)*sigma;
end

H_eq = sqrtm(K_n^(-1))*H*sqrtm(P);
switch Type
    case 'LMMSE'
     
        Phi = (H_eq'*H_eq+eye(N));
        G = Phi^(-1) *H'*K_n^(-1);
        SINR = 1./diag(Phi^(-1)) - 1;
        
    case 'MMSE_VBLAST'
        
        Phi = (H_eq'*H_eq+eye(N));
        G = Phi^(-1) *H'*K_n^(-1);

        SINR = zeros(N,1);
        index = [1:N]';
        for i = 1:N
            currentPhi = (H_eq'*H_eq+eye(size(index,1)));

            currentSINR = 1./diag(currentPhi^(-1))-1;
            [maxSINR, maxSINR_index] = max(currentSINR);
            
            SINR_index = index(maxSINR_index);
            SINR(SINR_index) = maxSINR;
           
            index(maxSINR_index)=[];
            H_eq(:,maxSINR_index)=[];
        end
        
end
end

