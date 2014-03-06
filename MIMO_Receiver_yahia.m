function [ x_rec,Phi ] = MIMO_Receiver_yahia( y,H,P,sigma, Type )
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
        
        Phi = (P*H'*K_n^(-1)*H+eye(N))^(-1);
        G = Phi *H'*K_n^(-1);
        
        x_rec = G*y;
        
    case 'MMSE_VBLAST'
        if size(sigma) == [M,M]
            K_n = sigma;
        elseif size(sigma,1) == M
            K_n = diag(sigma);
        else
            K_n = eye(M)*sigma;
        end
        
        
        
        H1=sqrtm(inv(K_n))*H;
        Tx_ant=size(H1,2);
        Rx=size(H1,1);
        current_H = H1;
        %current_y = y;
        %x = zeros(Tx_ant, 1);
        SNRlevels = zeros(Tx_ant, 1);
        x_index = [1:Tx_ant];
        for i_Tx = 1:Tx_ant
            MMSE_inv = inv(current_H'*current_H + 1 * eye(Tx_ant-i_Tx+1));
            [min_power, min_index] = min(diag(MMSE_inv));
            %             SNRlevels(x_index(min_index), 1) = 1/min_power-1;
            Phiindex(x_index(min_index), 1) = min_power;
            current_H = [current_H(:,1:min_index-1) current_H(:,min_index+1:(Tx_ant+1-i_Tx))];
            x_index = [x_index(1:min_index-1) x_index(min_index+1:(Tx_ant+1-i_Tx))];
            
        end
        %   x;
        x_rec=x_index;
        Phi=diag(Phiindex);
        
end
end

