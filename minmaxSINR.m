function [ P_op, gradient ] = minmaxSINR( H,P,sigma )

pNorm = 1; % 1: maximize sumRate
iterations = 50000;
seed_start=1000;
randn('state',seed_start);
rand('state',seed_start);
Px=rand(1,4);
Pi=Px.^2/(sum(Px.^2));
Pi1=4/(sum(Pi))*Pi;
P=diag(Pi1);
%P=diag([0.25 0.5 0.4 4-1.15]);

P=6*P;
[M,N] = size(H);
maxP = trace(P);
H_eq = sigma^(-1/2)*H;

Phi = P^(1/2)'*(H_eq'*H_eq)*P^(1/2)+eye(N);

gradient = zeros(N,iterations);
GradNorm = zeros(1,iterations);
sumRate = zeros(1,iterations);

X = sqrt(P);

SINRtgt =3*ones(N,1) ;%1./diag(Phi^(-1))-1;
% numerical Gradient
for j=1:iterations
    %iterate
    Phi = X'*(H_eq'*H_eq)*X+eye(N);
    Rate = real(log2(1./diag(Phi^(-1))));
    sumRate(j) = norm(Rate,pNorm);
    diffToTgt(j) = norm(1./diag(Phi^(-1))-1-SINRtgt);
    
    %calculate Gradient
    e = 10^(-4);
    diffToTgt_e = zeros(N,1);
    
    for i=1:N
        X_e = X; X_e(i,i) = X(i,i)+e;
%         X_e = X_e*sqrt(maxP)/sqrt(trace(X_e^2)); %normierung
        
        Phi_e = X_e'*(H_eq'*H_eq)*X_e+eye(N);
        diffToTgt_e(i) = norm(1./diag(Phi_e^(-1)) -1 -SINRtgt);
        
    end

    gradient(:,j) = (diffToTgt_e-diffToTgt(j))./e;
    
    step = .0001;
    X = X-step*diag(gradient(:,j));
%     if trace(X)>sqrt(maxP)
%         X = X*sqrt(maxP)/sqrt(trace(X^2));
%     end
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
    plot(diffToTgt)
end

Phi = P_op^(1/2)'*(H_eq'*H_eq)*P_op^(1/2)+eye(N);
end

