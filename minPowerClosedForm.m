function [ P_op ] = minPowerClosedForm( H,P,sigma )

[M,N,K] = size(H);
H_new = zeros(M,K,N);
for n=1:N
    for k=1:K
        H_new(:,k,n) = H(:,n,k);
    end
end
H = H_new;[M,N,K] = size(H);

Gamma(:,1) = ones(K,1); %target SINR
Gamma = diag(Gamma);

Omega_2 = @(k,j) (H(:,:,k)'*H(:,:,k))^(-1);
Omega_1 = @(k,j) Omega_2(k,j) *H(:,:,k)'*H(:,:,j)*H(:,:,j)'*H(:,:,k) * Omega_2(k,j);

F = zeros(K); noise = zeros(K,1);
for k=1:K
    for j=1:K
        if j~=k
            F(k,j) = max(eig(real(Omega_1(k,j))));
        else
            noise(k) = N*sigma*max(eig(Omega_2(k,j)));
        end
    end
end

P_star = real((eye(K)-Gamma*F)^(-1)*Gamma*noise);
P_op = diag(P_star);

end

