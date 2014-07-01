clear;

M = 2;
N = 1;
K = 2;

P = 10;

%Choose channel for iteration
%UPLINK CHANNEL!!!!!!
%H(:,:,k) represents uplink channel of kth user
H = 1/sqrt(2) * (randn(M,N,K) + 1i * randn(M,N,K));

%calculate DPC sum rate with iterative waterfilling
[capacity Covar] = iterative_waterfill(H,P,50)

%compute downlink covariance matrices
Downlink = MAC_to_BC(H, Covar)

