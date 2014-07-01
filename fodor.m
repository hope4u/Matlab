function [ P_op ] = fodor( H,P,sigma )
%FODOR implementation of the Algorithm for power Optimization with fairnes
%constraints form Fodors Paper
[M,N] = size(H);

K=1;% one antenna per User

%                 alpha = sqrt(P./K);
%                 chi = 1;
%                 d = 1;
% 
%                 R_n = 1/sigma * eye(N);
%                 R_z = 0; %inter cell interference
%                 R_v = R_z + R_n;
%                 R_H = alpha^2 * H'*R_v^(-1)*H;
%                 % G = (T'*R_H*T+eye(N))^(-1) * 1/alpha * T'*H'*R_v^(-1);
%                 T=eye(N);
%                 Phi = T'*R_H*T+eye(N);
%                 E = Phi^(-1);
% 
%                 gamma = 1./diag(Phi^(-1)) - 1;
% 
%                 gamma_ = 1/(N*max(eig((H'*H)^(-1))));


H_eq = sigma^(-1/2)*H*sqrtm(P);
Phi = (H_eq'*H_eq+eye(N));
E = Phi^(-1);


%% optimizing sum Power under SINR target
Gamma = .2*eye(N); % SINR Target

noise = sigma*1./diag(H'*H);  %%TODO: is this resonable since the result in the end is not SINR = Gamma but something lower...

p_star = diag(eye(N)*Gamma*noise);
Phi = []; E = []; noise = [];
for n = 1:N
    H_eq = sigma^(-1/2)*H(n,n)*sqrtm(P(n,n));
    Phi(n) = H_eq'*H_eq + eye(K);
    E(n) = Phi(n)^(-1);
    gamma(n) = 1/E(n) - 1;
    noise(n) = sigma*eig((H(n,n)'*H(n,n))^(-1));
end

p_star = eye(N)*Gamma*noise';
%% Optimal SINR traget selection



end

