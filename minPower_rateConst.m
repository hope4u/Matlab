function [ P_op, gradient ] = minPower_rateConst( H,P,sigma )

pNorm = 1; % 1: maximize sumRate
iterations = 1000;

[M,N] = size(H);
maxP = trace(P);
H_eq = sigma^(-1/2)*H;

Phi = P^(1/2)'*(H_eq'*H_eq)*P^(1/2)+eye(N);

gradient = zeros(N,iterations);
GradNorm = zeros(1,iterations);
sumRate = zeros(1,iterations);

X = sqrt(P);

minRate = 1;
% numerical Gradient
for j=1:iterations
    %iterate
    Phi = X'*(H_eq'*H_eq)*X+eye(N);
    Rate = real(log2(1./diag(Phi^(-1))));
%     sumRate(j) = norm(Rate,pNorm);
    Phiinv = Phi^(-1);
    for k=1:N
        sumRate(j) = sumRate(j) + (norm(Rate,pNorm)/N - Rate(k))^2;
    end
    sumRate(j) = sumRate(j)^(1/N);
%     diffToTgt(j) = norm(1./diag(Phi^(-1))-1-SINRtgt);
    
    %calculate Gradient
    e = 10^(-4);
    sumRate_e = zeros(N,1);
    
    for i=1:N
        X_e = X; X_e(i,i) = X(i,i)+e;
        X_e = X_e*sqrt(maxP)/sqrt(trace(X_e^2)); %normierung
        
        Phi_e = X_e'*(H_eq'*H_eq)*X_e+eye(N);
        Phiinv_e = Phi_e^(-1);
        
        Rate_e = real(log2(1./diag(Phi_e^(-1))));
%         sumRate_e(i) = norm(real(log2(1./diag(Phi_e^(-1)))),pNorm);
        for k=1:N
            sumRate_e(i) = sumRate_e(i) + (norm(Rate_e,pNorm)/N - Rate_e(k))^2;
        end
        
        sumRate_e(i) = sumRate_e(i)^(1/N);
    end

    gradient(:,j) = (sumRate_e-sumRate(j))./e;
    
    step = .001;
    X = X-step*diag(gradient(:,j));
    if trace(X^2)>maxP
        X = X*sqrt(maxP)/sqrt(trace(X^2));
    end
    P = X^2;
    P_now(:,j) = diag(P);
    
    GradNorm(j) = norm(gradient(:,j));
end
    
P_op = P;

if sigma == 1
    %plot Gradient and sumRate
    figure
    plot(GradNorm)
    hold all
    plot(sumRate)
end

Phi = P_op^(1/2)'*(H_eq'*H_eq)*P_op^(1/2)+eye(N);
end

