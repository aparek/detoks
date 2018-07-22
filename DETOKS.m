function [x,s,f,c,spindleVec,binary,binaryK,cost] = DETOKS(y,fs,fc,lam1,lam2,lam3,c1,c2)
%function [x,s,f,c,spindleVec,binary,binaryK,cost] = detoks(y,fs,d,fc,lam1,lam2,lam3,Nit, windowLength, c1,c2)
%Spindle detection using sparse regularization
%
% input - 
%         y - EEG channel
%         fs - sampling frequency 
%         fc - cut-off frequency of the low pass filter
%         lam1,lam2,lam3 - regularization parameters
%         c1,c2 - threshold for usage after applying TEAGER operator
%         
% output - 
%         x - transient component
%         s - oscillatory component
%         f - low-frequency component
%         c - STFT coefficients
%         spVec - bandpass filter applied to s
%         binS - binary vector containing the spindle location (1 - spindle, 0- not
%                   a spindle)
%         binK - binary vector containing the K-complex location
%         cost - cost function history for detoks
%
%
% Ankit Parekh
% Icahn School of Medicine at Mount Sinai
% (ankit.parekh@mssm.edu)
% Last EDIT: 7/21/2018
%
% Please cite as:  A. Parekh, I.W. Selesnick, D.M.Rapaport and I.Ayappa,
% Detection of K-complexes and sleep spindles (DETOKS) using sparse
% optimization, Journal of Neuroscience Methods, 251:37-46, Aug. 2015. 

mu = 0.5;                               % Scaled Augmented Lagrangian Parameter
d = 2;                                  % Degree of low pass filter
Nit = 30;                               % Number of iterations for Transient sep. algorithm
windowLength = fs;                      % DFT window length
y = y(:);                               % Convert to column vector
cost = zeros(Nit,1);                    % Cost function history
N = length(y);   
[A, B] = ABfilt(d, fc, N);              % Create the low-pass filter
H = @(x) A\(B*x);                       % H: high-pass filter
G = mu*(A*A') + 2*(B*B');               % G: Banded system
bn = nan + zeros(d, 1);     	        % bn : nan's to extend f to length N

% For Sleep spindle use [A1,A1H,~] = MakeTransforms('STFT',N, [1*fs 4 2 1*fs]);
[A1,A1H,~] = MakeTransforms('STFT',N, [2^nextpow2(windowLength) 4 2 2^nextpow2(windowLength)]);

c = A1H(y);
x = zeros(N,1);
d1 = zeros(N,1);
d2 = A1H(y);

b = (1/mu) * B' * ((A*A')\(B*y));
Ab = A1H(b);
wbar = waitbar(0, 'Running DETOKS');
for i = 1:Nit
    g1 = b + (x+d1);
    g2 = Ab + (c+d2);
    u1 = g1 - B' * (G\(B*(g1 + A1(g2)')));
    u2 = g2 - A1H(B' * (G\(B*(g1 + A1(g2)'))));
    x = soft(tvd(u1-d1,N,lam2/(mu)),lam1/(mu))';
    c = soft(u2-d2,lam3/(mu));
    d1 = d1 - (u1 - x);
    d2 = d2 - (u2-c);
    
    cost(i) = 0.5*sum(abs(H(y-x-real(A1(c)'))).^2) + ...
              lam1 * sum(abs(x)) + lam2 * sum(abs(diff(x))) + ... 
              lam3 * sum(abs(c(:))) ;
    waitbar(i/Nit,wbar, 'Running DETOKS...');
end
close(wbar);

s = real(A1(c)');
f = y - x - s - [bn; H(y-x-s); bn];       
[B0,A0] = butter(4, [11.5/(fs/2) 15.5/(fs/2)]);
spindleVec = filtfilt(B0,A0,s); 

binary = T(spindleVec)>c1;
f(isnan(f)) = 0;   

binaryK = T(f)>c2;

% The code snippet below is adapted from S.C. Warby et al. Nature Methods 
% 11, 385-392 (2014). 
E = binary(2:end)-binary(1:end-1);
sise = size(binary);

begins = find(E==1)+1;

if binary(1) == 1
    if sise(1) > 1
        begins = [1; begins];
    elseif sise(2) > 1
        begins = [1 begins];
    else
        error('The input signal is not one dimensional')
    end
elseif numel(begins) == 0 && binary(1) == 0
    begins = NaN;
end

ends = find(E==-1);
if binary(end) == 1
    if sise(1) > 1
        ends = [ends; length(binary)];
    elseif sise(2) > 1
        ends = [ends length(binary)];
    else
        error('The input signal is not one dimensional')
    end
elseif numel(ends) == 0 && binary(end) == 0
    ends = NaN;
end

Ek = binaryK(2:end)-binaryK(1:end-1);
sisek = size(binaryK);
beginsK = find(Ek==1)+1;

if binaryK(1) == 1
    if sisek(1) > 1
        beginsK = [1; beginsK];
    elseif sisek(2) > 1
        beginsK = [1 beginsK];
    else
        error('The input signal is not one dimensional')
    end
elseif numel(beginsK) == 0 && binaryK(1) == 0
    beginsK = NaN;
end

endsK = find(Ek==-1);
if binaryK(end) == 1
    if sisek(1) > 1
        endsK = [endsK; length(binaryK)];
    elseif sisek(2) > 1
        endsK = [endsK length(binaryK)];
    else
        error('The input signal is not one dimensional')
    end
elseif numel(endsK) == 0 && binary(end) == 0
    endsK = NaN;
end

[binary,~,~] = minimum_duration(binary,begins,ends,0.5,fs);
[binary,~,~] = maximum_duration(binary,begins,ends,3,fs);


[binaryK,~,~] = minimum_duration(binaryK,beginsK,endsK,0.3,fs);
[binaryK,~,~] = maximum_duration(binaryK,beginsK,endsK,3,fs);


%% Functions
% These functions are from S.C. Warby et al. Nature Methods 11, 385-392
% (2014)
function [DD,begins,ends] = minimum_duration(DD,begins,ends,min_dur,fs)
        % MINIMUM_DURATION - checks the sample duration of the spindles.
        % Input is a vector containing ones in the interval where the spindle is
        % and indexs describing the start and end of the spindle. The last two
        % inputs are the minimum duration given in seconds and the sampling
        % frequency given in Hz.
        % Output is a vector containing ones in the interval where the spindle with
        % duration longer than or equal to the minimum duration is and indexs
        % describing the start and end of the spindle.
        
        duration_samples = ends-begins+1;
        for k = 1:length(begins)
            if duration_samples(k) < min_dur*fs
                DD(begins(k):ends(k)) = 0;
                begins(k) = 0;
                ends(k) = 0;
            end
        end
        begins = begins(begins~=0);
        ends = ends(ends~=0);
end

function [DD,begins,ends] = maximum_duration(DD,begins,ends,max_dur,fs)
        % MAXIMUM_DURATION - checks the sample duration of the spindles.
        % Input is a vector containing ones in the interval where the spindle is
        % and indexs describing the start and end of the spindle. The last two
        % inputs are the maximum duration given in seconds and the sampling
        % frequency given in Hz.
        % Output is a vector containing ones in the interval where the spindle with
        % duration shorter than or equal to the maximum duration is and indexs
        % describing the start and end of the spindle.
        
        duration_samples = ends-begins+1;
        for k = 1:length(begins)
            if duration_samples(k) > max_dur*fs
                DD(begins(k):ends(k)) = 0;
                begins(k) = 0;
                ends(k) = 0;
            end
        end
        begins = begins(begins~=0);
        ends = ends(ends~=0);
end
end

