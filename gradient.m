function [ P_op ] = gradient( H,P,sigma )

[M,N] = size(H);
maxP = trace(P);
H_eq = sigma^(-1/2)*H;

Phi = P^(1/2)'*(H_eq'*H_eq)*P^(1/2)+eye(N);
Rate1 = real(log2(1./diag(Phi^(-1))))

%% analytical gradient

for j = 1:1000
    X = sqrt(P); X = X*sqrt(maxP)/sqrt(trace(X^2));
    Phi = X'*(H_eq'*H_eq)*X+eye(N);
    
    Gradient = zeros(N,1);
    for i = 1:N
        J=zeros(N); J(i,i)=1;
        Gradient(i) = real(log(2)*sum(diag(1./Phi^(-1)).*diag(Phi^(-1)*(J'*H'*1/sigma*H*X+X*H'*1/sigma*H*J)*Phi^(-1))));
    end
    
    X = X + .1*diag(Gradient);
    X = X*sqrt(maxP)/sqrt(trace(X^2));
    P = X^2;    

    for i=1:N
        plot(j,P(i,i))
        hold on
    end
end

%% numerical Gradient

% for j=1:1000
%     %iterate
%     X = sqrt(P); X = X*sqrt(maxP)/sqrt(trace(X^2));
%     Phi = X'*(H_eq'*H_eq)*X+eye(N);
%     Rate = real(log2(1./diag(Phi^(-1))));
% 
%     %calculate Gradient
%     e = 10^(-11);
%     Rate_e = Rate;
%     for i=1:N
%         X_e = X; X_e(i,i) = X(i,i)+e;
%         X_e = X_e*sqrt(maxP)/sqrt(trace(X_e^2)); %normierung
% 
%         Phi_e = X_e'*(H_eq'*H_eq)*X_e+eye(N);
%         Rate_e_i = real(log2(1./diag(Phi_e^(-1))));
%         Rate_e(i) = Rate_e_i(i);
%     end
% 
%     Gradient = (Rate_e-Rate)/e;
%     X = X+.1*diag(Gradient);
%     X = X*sqrt(maxP)/sqrt(trace(X^2));
%     P = X^2;
%     
%     for i=1:N
%         plot(j,P(i,i))
%         hold on
%     end
% end
    
P_op = P;

Phi = P_op^(1/2)'*(H_eq'*H_eq)*P_op^(1/2)+eye(N);
Rate2 = real(log2(1./diag(Phi^(-1))))
end

