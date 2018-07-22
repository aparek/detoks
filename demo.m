%% Spindle detection DEMO
% 
% This script provides a demo of the spindle detection used in [1]. The
% detection requires the use of several parameters, and at each step their
% significance is noted in this script. 
%
% Please cite as: [1] A. Parekh, I.W. Selesnick, D.M.Rapaport and I.Ayappa,
% Detection of K-complexes and sleep spindles (DETOKS) using sparse
% optimization, Journal of Neuroscience Methods, 251:37-46, Aug. 2015. 

%% Initialize

clear; close all; clc
warning('off','all');
format long g;
%% Load the EEG data
% EDF files
params.filename = 'excerpt2';
url = ['http://www.tcts.fpms.ac.be/~devuyst/Databases/DatabaseSpindles/',params.filename,'.edf'];
websave([params.filename,'.edf'],url);
[data, header] = lab_read_edf([params.filename,'.edf']);
N = length(data);
fs = header.samplingrate;   

% Visual detections
url = ['http://www.tcts.fpms.ac.be/~devuyst/Databases/DatabaseSpindles/Visual_scoring1_',params.filename,'.txt'];
websave(['Visual_scoring1_',params.filename,'.txt'],url);

url = ['http://www.tcts.fpms.ac.be/~devuyst/Databases/DatabaseSpindles/Visual_scoring2_',params.filename,'.txt'];
websave(['Visual_scoring2_',params.filename,'.txt'],url);

% Load the visual detections
fid = fopen(['Visual_scoring1_',params.filename,'.txt'],'r');
visualScorer1 = textscan(fid,'%f %f','headerlines',1,'Delimiter','%n');
visualScorer1 = [visualScorer1{1}, visualScorer1{2}];
fclose(fid);
vd1 = obtainVisualRecord(visualScorer1,fs,N);   

fid = fopen(['Visual_scoring2_',params.filename,'.txt'],'r');
visualScorer2 = textscan(fid,'%f %f','headerlines',1,'Delimiter','%n');
visualScorer2 = [visualScorer2{1}, visualScorer2{2}];
fclose(fid);
vd2 = obtainVisualRecord(visualScorer2,fs,N);  


%% Run the DETOKS algorithm

% Parameters
Hz = 1;                                                                     % Cut-off frequency for low-pass filter                                     
d = 2;                                                                      % Degree of low-pass filter
fc = Hz/(fs/2);                                                             % Normalize the cut-off frequency
c1 = 0.03; c2 = 1;                                                          % Cut-off value to be used after applying Teager operator. c1 - spindle, c2 - K
lam0 = 0.6;                                                                 % sparsity of the transients (higher value forces zero baseline)
lam1 = 7;                                                                   % sparsity of derivative of transients (higher value forces fewer steps)
lam2 = 8;                                                                   % sparsity of the oscillatory component (higher value forces zero baseline)
channel = 3;                                                                % CZ channel from excerpt 2
[x,s,f,c,spVec,binS,binK,cost] = DETOKS(data(channel,:),fs,fc,lam0,lam1,... % Run the DETOKS algorithm
                                            lam2,c1,c2);
% x - transient component
% s - oscillatory component
% f - low-frequency component
% c - STFT coefficients
% spVec - bandpass filter applied to s
% binS - binary vector containing the spindle location (1 - spindle, 0- not
% a spindle)
% binK - binary vector containing the K-complex location
% cost - cost function history for detoks

figure(1), clf
plot(cost, 'o-k')
grid on
box off
xlabel('Iteration')
ylabel('DETOKS Cost function')
                                        
%% Scoring

% In case scoring is needed: the following calculates the F1 score and
% other parameters based on the gold standard as visualScorer1 and
% visualScorer2
SCORE{2} = F1score(binS,vd1,vd2);                                           % Calculate the scoring 
SCORE{2}{1}'                                                                % Display the score type
SCORE{2}{2}                                                                 % Display the values
%% Plot the result

epoch = 1;                                                                  % Choose an epoch to view
n = 0:30*fs-1;
gap = 125;
fig1 = figure(2); clf;
set (fig1, 'Units', 'normalized', 'Position', [0,0,1,1]);
plot(n/fs, data(channel, (epoch-1)*30*fs+1:epoch*30*fs), ...
     n/fs, x((epoch-1)*30*fs+1:epoch*30*fs) - gap,...
     n/fs, f((epoch-1)*30*fs+1:epoch*30*fs) - 2*gap,...
     n/fs, s((epoch-1)*30*fs+1:epoch*30*fs)*5 - 3*gap,...
     n/fs, binS((epoch-1)*30*fs+1:epoch*30*fs)*50 - 4*gap);
box off
grid on
ylim([-5*gap 2*gap])
set(gca,'YTick',[])
legend('Raw EEG','Transient','Low-frequency','Oscillatory (x5)','Spindle')