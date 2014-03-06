function [ y,H ] = MIMO_Channel( x, M,N,sigma )

H = 1/sqrt(2)*(randn(M,N)+1i*randn(M,N));
v = sigma*1/sqrt(2)*(randn(M,1)+1i*randn(M,1));

y = H*x+v;

end

