%%  This file provides a demo of the DETOKS algorithm
%   
%   Please cite as: 
%   A. Parekh, I.W. Selesnick, D. Rapaport, I.Ayappa, Detection of
%   K-complexes and sleep spindles (DETOKS) using sparse optimization,
%   Journal of Neuroscience Methods, vol. 251, pp. 37-46, 2015. 
%   
%   Last Edit: March 29th, 2016. Ankit Parekh(ankit.parekh@nyu.edu)
%   Copyright (c) 2015. 
%% Initialization

clear; close all; clc;
%% Create synthetic signal

% Initialize noise level and signal length
rng('default')
sigma = 0.1;
fs = 25;
N = 10*fs;
n = 0:N-1;


% Generate synthetic components
% low-frequency
f = 0.25;
f = sin(2*f/fs*pi*n);      

% oscillatory
s = zeros(size(n));
s(100+(1:30)) = sin(2*pi*0.15*(1:30)) .* hamming(30)';
x = zeros(size(n));

% sparse piece-wise constant
x(50:60) = -1;
x(140:160) = 1;

w = sigma * randn(size(n));
y = f+s+x + w;

figure(1), clf
gap = 3;
plot(n/fs, f+s+x, n/fs, x-gap, n/fs, f-2*gap, n/fs, s-3*gap, n/fs, w-4*gap)
box off
grid on
legend('Clean Signal','Transient Component', 'Low-frequency Component', 'Oscillatory Component','AWGN','location','northoutside')
ylim([-4.5*gap gap])
set(gca,'YTick',...
    [-3*gap -2*gap  -gap 0],'YTickLabel',0)
xlabel('Time (s)')

%% Run DETOKS to estimate the components

% Parameters
Hz = 0.5;
lam1 = 0.01;
lam2 = 0.3;
lam3 = 0.15;
d = 2;      
fc = Hz/(fs/2);
Nit = 30;
mu = 0.8;

time = cputime;
[x,s,f,cost] = DETOKS(y,fs,d,fc,lam1,lam2,lam3,Nit,mu);
time = cputime-time;

figure(2), clf
plot(cost, 'k')
box off
xlabel('Iteration')
ylabel('Cost function value')
title(sprintf('Cost function history. Time taken = %1.2f seconds', time))

%% Plot the components
figure(3), clf

residual = y-(x+f+s)';
gap = 3;
plot(n/fs, y, n/fs, x-gap, n/fs, f-2*gap, n/fs, s-3*gap, n/fs, residual-4*gap)
box off
grid on
legend('Noisy Input Signal','Transient Component', 'Low-frequency Component', 'Oscillatory Component','Residual','location','northoutside')
title('Estimated components using DETOKS')
ylim([-4.5*gap gap])
set(gca,'YTick',...
    [-3*gap -2*gap  -gap 0],'YTickLabel',0)
xlabel('Time (s)')
