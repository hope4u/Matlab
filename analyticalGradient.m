function [ P_op,gradient ] = analyticalGradient( H,P,sigma )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


pNorm = 1; % 1: maximize sumRate
iterations = 1000;

[M,N] = size(H);
maxP = trace(P);
H_eq = sigma^(-1/2)*H;

Phi = P^(1/2)'*(H_eq'*H_eq)*P^(1/2)+eye(N);

sumRate = zeros(1,iterations);
gradient = zeros(N,iterations);
GradNorm = zeros(1,iterations);


% ableitung der Norm
%iterate
X = sqrt(P);

for n=1:iterations
    
    Phi = X'*(H_eq'*H_eq)*X+eye(N);
    invPhi = Phi^(-1);
    Rate = real(-log2(diag(Phi^(-1))));
    sumRate(n) = norm(Rate,pNorm);
    
    for j=1:N
        E = zeros(N); E(j,j) = 1;
        DRate = zeros(N,1);
        
            
        Dtrace = trace(2*E*X);
        Dxhhx = X*H'*H*E + E*H'*H*X;
        DPhi = maxP/(trace(X^2))^2 * (Dxhhx*trace(X^2) - X*H'*H*X*Dtrace);

        DinvPhi = -Phi^(-1) * DPhi * Phi^(-1);
            
        for i=1:N
            DRate(i) = real(-1/(log(2)*invPhi(i,i)) * DinvPhi(i,i));
        end
        
        gradient(j,n) = (sum(Rate.^pNorm)^(1/pNorm-1)) * sum(Rate.^(pNorm-1) .* DRate);
    
    end
    
    X = X +.1*diag(gradient(:,n));
    X = X*sqrt(maxP)/sqrt(trace(X^2));

    P = X^2;    

    GradNorm(n) = norm(gradient(:,n));
end

P_op = P;

if sigma == 1
    %plot Gradient and sumRate
    figure
    plot(GradNorm)
    hold on
    plot(sumRate)
end

Phi = P_op^(1/2)'*(H_eq'*H_eq)*P_op^(1/2)+eye(N);
end

