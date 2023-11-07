function Out_MNSync = MNSync_ZHY(y,m,mu,Mpr,re)
%
% This function calculates multipleresolution nonsynchronous entropy (MNSyncEn) of a synchronous series y
%
% Inputs:
% y: synchronous series - a vector of size 1 x N (the number of sample points)
% m: embedding dimension
% mu: number of intervals
% Mpr: number of samples per revolution
% re: resolution coefficient

% Outputs:
% Out_MNSyncEn: MNSyncEn of the synchronous series y (complexity of the nonsynchronous component in the synchronous series y) 

% Please cite the paper entitled "Multiresolution nonsynchronous entropy:
% Measurement approach for synchronous series analysis and its application
% in fault diagnosis of rotating machinery," if you use this code.
%

% Yanqing ZHAO
% zhaoyanqing@hyit.edu.cn
%2022.10.21
% HYIT
%%

Out_MNSync = zeros(1,re);
for ii =1:re
P =zeros(1,mu);

NN = length(y);
M= NN-(m-1)*ii*Mpr;  % Number of vector sequences in the phase space
X_1 = zeros(M,m);
X = zeros(M,m);
for i=1:M      % Phase space reconstruction and Mutipleresolutions zero-centering techneque   
    for j=1:m
        X(i,j)=y(i+(j-1)*ii*Mpr);
    end
    X_1(i,:) = X(i,:)-1/m*sum(X(i,:));       
end

D_D = zeros(1,M-1);
A = zeros(1,M-1);
B = zeros(1,M-1);
C = zeros(1,M-1);
for i=1:M-1 

     A(i) = sum(X_1(i,:).*X_1(i+1,:));
     B(i) = (sum(X_1(i,:).^2)).^0.5;
     C(i) = (sum(X_1(i+1,:).^2)).^0.5;
     D_D(i) =A(i)/(B(i)*C(i)); 
end

Partition = -1+1/mu:2/mu:1;
N=hist(D_D,Partition);
N_sum = sum(N);

P = N./N_sum;
P_1 = P(find(P~=0));
Out_MNSync(ii) = -1/log(mu)*sum(P_1.*log(P_1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end