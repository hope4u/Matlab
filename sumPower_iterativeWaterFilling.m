function [ H_op,P_op ] = sumPower_iterativeWaterFilling( H,P,sigma )
% initialize
M = size(H,1);
K = size(H,2);
N=1; % umber of Antennas per User

%% Algorithm 1
% n=K-1;
% Q = zeros(K,n);
% Q(:,[-(K-2)+n:n]) = trace(P)/(K*N);
% 
% for n=n+1:n+10;
%     for i=1:K
%         Z = eye(M);
%         for j=1:K-1
%             ind = mod(i+j-1,K)+1;
%             Z = Z + H(:,ind)*Q(ind,n-K+j)*H(:,ind)';
%         end
%         G(:,i,n) = H(:,i)'*Z^(-1/2);
% 
%         [U(i),D(i),V(i)]=svd(G(:,i,n)'*G(:,i,n));
%     end
% 
%     u=(trace(P)+sum(D.^-1))/K;
%     A = u-D.^(-1);
% 
%     Q(:,n) = U.*A.*conj(U);
% end

%% Algorithm 2
n=1;
Q = zeros(K,n);
Q(:,n) = trace(P)/(K*N);

for n=n+1:n+50;
    for i=1:K
        Z = eye(M);
        for j=1:K-1
            ind = mod(i+j-1,K)+1;
            Z = Z + H(:,ind)*Q(ind,n-1)*H(:,ind)';
        end
        G(:,i,n) = H(:,i)'*Z^(-1/2);

        [U(i),D(i),V(i)]=svd(G(:,i,n)'*G(:,i,n));
    end

    u=(trace(P)+sum(D.^-1))/K;
    A = u-D.^(-1);
    S(:,n) = U.*A.*conj(U);
    
    Q(:,n) = 1/K * S(:,n)+(K-1)/K*Q(:,n-1);
end

%% return
P_op = diag(Q(:,n));
H_op = H;

end

