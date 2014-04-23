function [ P_op ] = fodorPrecodingOptimization( H,P,sigma )
[M,N] = size(H);

t=0; K=1;
epsilon(1) = 1;
P_tot = trace(P);
Gamma(:,1) = ones(N,1);

for k=1:M
    T(:,k) = ones(N,1);
end

for t=2:100
    for k = 1:K

        c(:,k) = diag( ...
                (H'*(N*sigma*eye(M))*H+1/P_tot*eye(N))^(-1));
        
        T_new(:,k) = sqrt(c(:,k)*N./sum(c(:,k)));
        
        Gamma(:,t) = Gamma(:,t-1).*epsilon(t-1);
        P(:,k,t) = c(:,k)./abs(T_new(:,k)).^2 .*(Gamma(k,t) + 1);
        P_tot = sum(P(:,k,t));
    end
    
    epsilon(t) = epsilon(t-1);
end

P_op = diag(T_new.*P(:,k,t));

end

