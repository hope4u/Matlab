%% calculate SINR

P_op
H_eq = sigma^(-1/2)*H*sqrtm(P_op)
Phi = (H_eq'*H_eq+eye(N))
SINR = 1./diag(Phi^(-1)) - 1

%% calculate Rate

P_op
H_eq = sigma^(-1/2)*H*sqrtm(P_op)
Phi = (H_eq'*H_eq+eye(N))
Rate = log2(1./diag(Phi^(-1)))

%% calculate sumRate

P_op
H_eq = sigma^(-1/2)*H*sqrtm(P_op)
Phi = (H_eq'*H_eq+eye(N))
Rate = sum(log2(1./diag(Phi^(-1))))