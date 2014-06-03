function [ P_op, gradient, diffToTgt ] = powermin( H,P,sigma )

iterations = 200;


[M,N] = size(H);
H_eq = sigma^(-1/2)*H;

Phi = P^(1/2)'*(H_eq'*H_eq)*P^(1/2)+eye(N);

gradient = zeros(N,iterations);
GradNorm = zeros(1,iterations);
diffToTgt = zeros(1,iterations);

X = sqrt(P);

SINRtgt =ones(N,1) ;%1./diag(Phi^(-1))-1;
step=.1;
% numerical Gradient
for j=1:iterations
    %iterate
    Phi = X'*(H_eq'*H_eq)*X+eye(N);
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
    
%     if trace(X)>sqrt(maxP)
%         X = X*sqrt(maxP)/sqrt(trace(X^2));
%     end


    step_old = step(end); step=step_old;
    diffToTgt_new=[]; fin = 0;
    for a=1:20
 
        X_new = X-step(a)*diag(gradient(:,j));

        if fin == 1
            X=X_new;
            break
        end
        % adaptive StepSize
        Phi = X_new'*(H_eq'*H_eq)*X_new+eye(N);
        diffToTgt_new(a) = norm(1./diag(Phi^(-1))-1-SINRtgt);

        if diffToTgt_new(a)<diffToTgt(j)
            if a>1 && diffToTgt_new(a)>diffToTgt_new(a-1)
                step(a+1)=step(a-1);
                fin = 1;
            else
                if a==20
                    X=X_new;
                end
                step(a+1)=step(a)*10;
            end
        else
            if a>1 && diffToTgt(j)>diffToTgt_new(a-1)
                step(a+1) = step(a-1);
                fin = 1;
            else
                step(a+1) = step(a)*1/10;
            end
        end
    end

    P = X^2;
    P_now(:,j) = diag(P);    
    GradNorm(j) = norm(gradient(:,j));
    if diffToTgt(j)<.001
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
%     plot(diffToTgt)
% end

Phi = P_op^(1/2)'*(H_eq'*H_eq)*P_op^(1/2)+eye(N);
end

