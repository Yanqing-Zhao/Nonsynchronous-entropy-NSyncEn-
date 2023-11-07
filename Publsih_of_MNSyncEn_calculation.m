% The Matlab code is used to calculat MNSyncEn values. We carried out ten 
% trials to reduce the randomness for each synchronous series.

clear;
clc;

%%Stable
% load('Milling_stable_data.mat')
% x_singal = 1000*x_singal; % Unit of length (mm)

% Hopf bifurcation
% load('Milling_Hopf_bifurcation_data.mat')
% x_singal = 1000*x_singal; % Unit of length (mm)

% Period-2 bifurcation
% load('Milling_period-2_bifurcation_data.mat')
% x_singal = 1000*x_singal; % Unit of length (mm)

% Period-3 bifurcation
load('Miling_period-3_bifurcation_data.mat')
x_singal = 1000*x_singal; % Unit of length (mm)

N_N = length(x_singal);
Entropy = zeros(N_N,8);

j = 10;
m =4;          % Embedding dimension
mu = 30;       % Number of intervals
Mpr = 400;     % Number of samples per revolution
re = 3;        % Resolution coefficient
Out_MNSyncZHY=zeros(re,j);
Out_MNSyncLC=zeros(re,j);

for ii =1:j   % We carried out ten trials to reduce the randomness for each synchronous series.
 
yy2 = awgn(x_singal,5,'measured'); % Signal-to-noise ratio 5dB
yy  = zscore(yy2(:));

Out_MNSyncZHY(1:re,ii) = MNSync_ZHY(yy,m,mu,Mpr,re); 
display(ii)
end

En1 = Out_MNSyncZHY';
MNSyncZHY_mean = zeros(re,1);
MNSyncZHY_mean(1,:) = mean(En1(:,1));
MNSyncZHY_mean(2,:) = mean(En1(:,2));
MNSyncZHY_mean(3,:) = mean(En1(:,3));
display(Out_MNSyncZHY)
display(MNSyncZHY_mean);

