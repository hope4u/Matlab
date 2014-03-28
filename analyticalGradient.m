function [ P_op ] = analyticalGradient( H,P,sigma )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


pNorm = 1; % 1: maximize sumRate
iterations = 100;

[M,N] = size(H);
maxP = trace(P);
H_eq = sigma^(-1/2)*H;

Phi = P^(1/2)'*(H_eq'*H_eq)*P^(1/2)+eye(N);
Rate1 = real(log2(1./diag(Phi^(-1))))


Gradient = zeros(N,iterations);
GradNorm = zeros(1,iterations);
sumRate = zeros(1,literations);

% analytical Gradient
for j=1:iterations
    %iterate
    X = sqrt(P); X = X*sqrt(maxP)/sqrt(trace(X^2));
    Phi = X'*(H_eq'*H_eq)*X+eye(N);
    Rate = real(log2(1./diag(Phi^(-1))));
    sumRate(j) = norm(Rate,pNorm);
    
    %calculate Gradient
    e = 10^(-6);
    sumRate_e = sumRate(j);
    
    for i=1:N
        X_e = X; X_e(i,i) = X(i,i)+e;
        X_e = X_e*sqrt(maxP)/sqrt(trace(X_e^2)); %normierung

        Phi_e = X_e'*(H_eq'*H_eq)*X_e+eye(N);
        sumRate_e(i) = norm(real(log2(1./diag(Phi_e^(-1)))),pNorm);
        
    end

    Gradient(:,j) = (sumRate_e-sumRate(j))./e;
    X = X+.1*diag(Gradient(:,j));
    X = X*sqrt(maxP)/sqrt(trace(X^2));
    P = X^2;
    
    GradNorm(j) = norm(Gradient(:,j));
end
    
P_op = P;

if sigma == 1
    plot Gradient and sumRate
    figure
    plot(GradNorm)
    plot(sumRate)
end

Phi = P_op^(1/2)'*(H_eq'*H_eq)*P_op^(1/2)+eye(N);
Rate2 = real(log2(1./diag(Phi^(-1))))



end

