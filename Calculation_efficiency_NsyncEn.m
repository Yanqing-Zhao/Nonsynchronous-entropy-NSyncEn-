
clear
clc
close all


Mpr = 200;
K = 6;
m =4;
mu = 100;
j = 1000;
Out_Nsync=zeros(1,j);
%% calculation time is counted in computing 1000*5 entropy values 
% x1 = load('Delta1p3.mat');   % 1000 entropy values 
% x1 = load('Delta1p35.mat');  % 1000 entropy values 
% x1 = load('Delta1p4.mat');   % 1000 entropy values 
% x1 = load('Delta1p45.mat');  % 1000 entropy values 
x1 = load('Delta1p5.mat');     % 1000 entropy values 
x2 = getfield(x1,'xx');
y1 = x2(1:10:end); % Set the number of samples per excitation period to 200

%% Calculation cost test ---NSyncEn
tic
for ii =1:j         % generate 1000 NSyncEn

y2 = awgn(y1,15);   % SNR=15dB
y3  = zscore(y2(:));



Out_Nsync(ii) = NSyncEn(y2,Mpr,K,m,mu);


end
toc
