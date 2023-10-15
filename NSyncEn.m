function Out_NSyncEn = NSyncEn(y,Mpr,K,m,mu)
%
% This function calculates nonsynchronous entropy (NSyncEn) of a synchronous series y
%
% Inputs:
% y: synchronous series - a vector of size 1 x N (the synchronous series length)
% Mpr: number of samples per revolution
% K: number of revolutions
% m: embedding dimension
% mu: resolution coefficient


% Outputs:
% Out_NSyncEn: NSyncEn of the synchronous series y (complexity of the nonsynchronous component in the synchronous series y) 


% Yanqing ZHAO
% zhaoyanqing@hyit.edu.cn
%2022.10.21
% Huaiyin Institute of Technology (HYIT)
%%

P =zeros(1,mu);

M=(K-m)*Mpr;                       %Number of vector sequences in the phase space
X_1 = zeros(M,m);
X = zeros(M,m);
for i=1:M                           %Phase space reconstruction
    for j=1:m
        X(i,j)=y(i+(j-1)*Mpr);
    end
    X_1(i,:) = X(i,:)-1/m*sum(X(i,:));  %Zero-centering techneque   
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
Out_NSyncEn = -1/log(mu)*sum(P_1.*log(P_1));


end