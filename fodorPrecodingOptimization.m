function [ P_op ] = fodorPrecodingOptimization( H,P,sigma )

[M,N,K] = size(H);

t=0;
epsilon(1) = 1;
kappa=.1;
Gamma(:,1) = ones(K,1); %target SINR

P_tot = trace(P)/N; P = P_tot;

T = ones(N,K,1);

for t=2:1000
    for k = 1:K

        int = 0;
        for j=1:K
            if j~=k
                int = int + P(j,t-1)*H(:,:,j)*diag(T(:,j,t-1))*diag(T(:,j,t-1))'*H(:,:,j)';
            end
        end       
        
        c(:,k,t) = real(diag( ...
            (H(:,:,k)'*(int+N*sigma*eye(M))^(-1)*H(:,:,k)+1/P(k,t-1)*eye(N))^(-1)));
        
        T(:,k,t) = sqrt(N*c(:,k,t)./sum(c(:,k,t)));
        
        Gamma(k,t) = real(Gamma(k,1).*epsilon(t-1));
        P(k,t) = sum(c(:,k,t))/N*(Gamma(k,t) + 1);
    end
    
    epsilon(t) = epsilon(t-1); %power Optimization
%     epsilon(t) = max(0,epsilon(t-1)-kappa*(sum(P(:,t))-P_tot)); %throughput maximization
end

for k=1:K
    P_op(:,:,k) = diag(T(:,k,t))*P(k,t)*diag(T(:,k,t))';
end

end