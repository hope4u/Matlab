function [ P_op ] = fodorPrecodingOptimization2( H,P,sigma )
[M,N] = size(H);

t=0; K=N;
epsilon(1) = 1;
kappa=.001;
Gamma(:,1) = ones(N,1)*.1; %target SINR

P_tot = trace(P);
P = diag(P);

T = ones(K,1);

for t=2:1000
    for k = 1:K

        int = 0;
        for j=1:K
            if j~=k
                int = int + P(j,t-1)*H(:,j)*T(j)*T(j)'*H(:,j)';
            end
        end       
        
        
        c(k) = diag( ...
                (H(:,k)'*(int+sigma*eye(M))^(-1)*H(:,k)+1/P(k,t-1)*eye(1))^(-1));
        
        Gamma(:,t) = Gamma(:,t-1).*epsilon(t-1);
        P(k,t) = c(k)*(Gamma(k,t) + 1);
    end
    
    epsilon(t) = epsilon(t-1);
    epsilon(t) = epsilon(t-1)-kappa*(sum(P(:,t))-P_tot);
    if epsilon(t)<0
        epsilon(t) = 0;
    end
end

P_op = diag(P(:,t))

end