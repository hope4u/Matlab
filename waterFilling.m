function [ H_wf,P_wf,U,S,V ] = waterFilling(H,P,sigma)
N = size(H,2);
[U,S,V]=svd(H);
s = diag(S'*S);
v = trace(P);
for k=1:N
    if v > sigma/s(k)
        v = (sum(sigma./s(1:k)) + trace(P))/k;
    end
end
P_wf = v - sigma./(S'*S);
P_wf(P_wf<0)=0;
P_wf=V*P_wf*V';
H_wf=H;

end