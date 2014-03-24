function [ P_op ] = gradient( H,P,sigma )
[M,N] = size(H);

H_eq = sigma^(-1/2)*H;

    for i=1:50
    Phi = P^(1/2)'*(H_eq'*H_eq)*P^(1/2)+eye(N);

    constGrad = [ones(N,1) eye(N)];

    n = diag(H_eq'*H_eq);
    nx = diag((P^(1/2)'*(H_eq'*H_eq)*P^(1/2))^(-1));
    rateGrad = - n ./ (nx * log(2)+log(2));

    s=-(eye(N)-constGrad*(constGrad'*constGrad)^(-1)*constGrad')*rateGrad;
    alpha = - .1*sum(log2(1./diag(Phi^(-1)))) / (s'*rateGrad);

    P = real(diag(diag(P) + alpha*s));
    end
    P_op = P;
end

