function [ H_op,P_op ] = sumPower_iterativeWaterFilling( H,P,sigma )


%% initialize
M = size(H,1);
K = size(H,2);
N=1; % umber of Antennas per User

n=K-1;
Q = zeros(K,n);
Q(:,[-(K-2)+n:n]) = trace(P)/(K*N);

for n=n+1:n+10;
    for i=1:K
        Z = eye(M)*sigma;
        for j=1:K-1
            ind = mod(i+j-1,K)+1;
            Z = Z + H(:,ind)*Q(ind,n-K+j)*H(:,ind)';
        end
        G(:,i,n) = H(:,i)'*Z^(-1/2);

        [U(i),D(i),V(i)]=svd(G(:,i,n)'*G(:,i,n));
    end

    u=(trace(P)+sum(D.^-1))/K;
    A = u-D.^(-1);

    Q(:,n) = U.*A.*conj(U);
end
P_op = diag(Q(:,n));
trace(P_op)
H_op = H;

end

