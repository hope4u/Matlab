function [ P_op,gradient, diffToTgt ] = Copy_of_analyticalGradient( H,P,sigma )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


iterations = 5000000;

[M,N] = size(H);
H_eq = sigma^(-1/2)*H;
tgtSINR = .1*ones(N,1);

Phi = P^(1/2)'*(H_eq'*H_eq)*P^(1/2)+eye(N);

diffToTgt = zeros(1,iterations);
gradient = zeros(N,iterations);
GradNorm = zeros(1,iterations);


% ableitung der Norm
%iterate
X = sqrt(P);

for n=1:iterations
    
    Phi = X'*(H_eq'*H_eq)*X+eye(N);
    invPhi = Phi^(-1);

    diffToTgt(n) = norm((1./diag(Phi^(-1))-1-tgtSINR));
    
    for j=1:N
        E = zeros(N); E(j,j) = 1;
        
            
        Dxhhx = X*H'*H*E + E*H'*H*X;
%         DPhi = maxP/(trace(X^2))^2 * (Dxhhx*trace(X^2) - X*H'*H*X*Dtrace);
        DPhi = Dxhhx;
        DinvPhi = -Phi^(-1) * DPhi * Phi^(-1);
            
%         for i=1:N
% %             DRate(i) = real(-1/(log(2)*invPhi(i,i)) * DinvPhi(i,i));
% %             DSINRtot = 1/invPhi^2 * DinvPhi;
%             DSINR(i) = -(invPhi(i,i))^(-2) * DinvPhi(i,i);
%         end
        
%         gradient(j,n) = sum((Rate-tgtRate).^2)^(-1/2) * sum((Rate-tgtRate) .* DRate);
        gradient(j,n) = real(sum((1./diag(invPhi)-1-tgtSINR).^2)^(-1/2) * sum((1./diag(invPhi)-1-tgtSINR) .* -1./(diag(invPhi).^2) .* diag(DinvPhi)));
    end
    
    step = .000001;
    X = X-step*diag(gradient(:,j));
%     X = X*sqrt(maxP)/sqrt(trace(X^2));

    P = X^2;    

    GradNorm(n) = norm(gradient(:,n));
end

P_op = P;

if sigma == 1
    %plot Gradient and sumRate
    figure
    plot(GradNorm)
    hold all
    plot(diffToTgt)
end

Phi = P_op^(1/2)'*(H_eq'*H_eq)*P_op^(1/2)+eye(N);
end

