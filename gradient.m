function [ P_op ] = gradient( H,P,sigma )
[M,N] = size(H);

H_eq = sigma^(-1)*H*P^(1/2);

Phi = H_eq'*H_eq+eye(M);



end

