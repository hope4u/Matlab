function [ P_op ] = fodorloops( H,P, sigma )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[M,N] = size(H);
K = M;


Omega_2 = @(k,j) (H(k,k)'*H(k,k))^(-1);
Omega_1 = @(k,j) Omega_2(k,j) * H(k,k)'*H(k,j)*H(k,j)'*H(k,k) * Omega_2(k,j);

F = zeros(K,K);
noise = zeros(K,1);
for k=1:K
    R_z = 0;
    R_n = sigma;

    noise(k) = sigma * eig(Omega_2(k,k)); 
    for j=1:K
        if(k~=j)
            F(k,j) = eig(Omega_1(k,j));
            R_z = R_z + P(j,j)*H(k,j)*H(k,j)';
        end
    end
    R_v = R_z+R_n;
    R_H = P(k,k)*H(k,k)'*R_v^(-1)*H(k,k);
    E(k) = (1+R_H)^(-1);
    
    SINR = 1./E -1;
end

Gamma = eye(K)*(.2);

p_star = (eye(K)-Gamma*F)^(-1) *Gamma*noise;
P_op = diag(p_star);
end

