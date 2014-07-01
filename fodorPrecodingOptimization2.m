function [ P_op ] = fodorPrecodingOptimization2( H,P,sigma,TgtRate )

[M,N,K] = size(H);

H_new = zeros(M,K,N);
for n=1:N
    for k=1:K
        H_new(:,k,n) = H(:,n,k);
    end
end
H = H_new;[M,N,K] = size(H);

t=0;
epsilon(1) = 1;
kappa=.1;
P_tot = trace(P); P = diag(P);
Gamma(:,1) = ones(K,1); %target SINR


T = ones(N,K,1);

Gamma(:,1) = ones(K,1)*(2^TgtRate-1); %target SINR
SINRtgt = Gamma(:,1);



for t=2:10000
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
        
        SINR_now(k,t-1)= P(k,t-1)/c(:,k,t) -1;
    end
    
    epsilon(t) = epsilon(t-1); %power Optimization
%     epsilon(t) = max(0,epsilon(t-1)-kappa*(sum(P(:,t))-P_tot)); %throughput maximization
    if abs(P(k,t-1)-P(k,t))<.0001
        break;
    end
end

for k=1:K
    P_op(k,k,1) = diag(T(:,k,t))*P(k,t)*diag(T(:,k,t))';
end

t=t+1;
for k=1:K
        int = 0;
        for j=1:K
            if j~=k
                int = int + P(j,t-1)*H(:,:,j)*diag(T(:,j,t-1))*diag(T(:,j,t-1))'*H(:,:,j)';
            end
        end       
        
        c(:,k,t) = real(diag( ...
            (H(:,:,k)'*(int+N*sigma*eye(M))^(-1)*H(:,:,k)+1/P(k,t-1)*eye(N))^(-1)));
        
        SINR_now(k,t-1)= P(k,t-1)/c(:,k,t) -1;
        P_now(k,:)=P(k,:);
end  


end