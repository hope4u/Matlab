function [ P_op, gradient ] = minmaxSINR( H,P,sigma )

pNorm = -100; % 1: maximize sumRate
iterations = 100;

[M,N] = size(H);
maxP = trace(P);
H_eq = sigma^(-1/2)*H;

Phi = P^(1/2)'*(H_eq'*H_eq)*P^(1/2)+eye(N);

gradient = zeros(N,iterations);
GradNorm = zeros(1,iterations);
sumRate = zeros(1,iterations);

X = sqrt(P);
% numerical Gradient
step = .1;
for j=1:iterations
    %iterate
    Phi = X'*(H_eq'*H_eq)*X+eye(N);
    Rate = real(log2(1./diag(Phi^(-1))));
    sumRate(j) = norm(Rate,pNorm);
    
    %calculate Gradient
    e = 10^(-4);
    sumRate_e = zeros(N,1);
    
    for i=1:N
        X_e = X; X_e(i,i) = X(i,i)+e;
        X_e = X_e*sqrt(maxP)/sqrt(trace(X_e^2)); %normierung
        
        Phi_e = X_e'*(H_eq'*H_eq)*X_e+eye(N);
        sumRate_e(i) = norm(real(log2(1./diag(Phi_e^(-1)))),pNorm);
        
    end

    gradient(:,j) = (sumRate_e-sumRate(j))./e;

    %step_old = step(end); 
    step_old = step(1); step=step_old;
    sumRate_new=[]; fin = 0;
    for a=1:40
 
        X_new = X+step(a)*diag(gradient(:,j));
        X_new = X_new*sqrt(maxP)/sqrt(trace(X_new^2));

        if fin == 1
            X=X_new;
            break
        end
        % adaptive StepSize
        Phi = X_new'*(H_eq'*H_eq)*X_new+eye(N);
        Rate = real(log2(1./diag(Phi^(-1))));
        sumRate_new(a) = norm(Rate,pNorm);

        if sumRate_new(a)>sumRate(j)
            if a>1 && sumRate_new(a)<sumRate_new(a-1)
                step(a+1)=step(a-1);
                fin = 1;
            else
                if a==20
                    X=X_new;
                end
                step(a+1)=step(a)*10;
            end
        else
            if a>1 && sumRate(j)<sumRate_new(a-1)
                step(a+1)=step(a-1);
                fin = 1;
            else
                step(a+1) = step(a)*1/10;
            end
        end
    end
    P = X^2;
    GradNorm(j) = norm(gradient(:,j));
    if GradNorm(j) < .001
        break
    end
    if j>1 && ~sum(gradient(:,j-1)-gradient(:,j))
        break
    end
end


P_op = P;

% if sigma == 1
%     %plot Gradient and sumRate
%     figure
%     plot(GradNorm)
%     hold all
%     plot(sumRate)
% end

Phi = P_op^(1/2)'*(H_eq'*H_eq)*P_op^(1/2)+eye(N);
end

