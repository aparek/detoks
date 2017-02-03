function [A, B, B1, D, a, b, b1] = ABfilt(deg, fc, N)
% [A, B, B1] = ABfilt(d, fc, N)
%
% Banded matrices for zero-phase high-pass recursive filter.
% The matrices are created as 'sparse' structures.
%
% INPUT
%   d  : degree of filter is 2d
%   fc : cut-off frequency (normalized frequency, 0 < fc < 0.5)
%   N  : length of signal
%
% OUTPUT
%   A, B, B1 : banded filter matrices
%       with B = B1*D where D is the first-order difference matrix
%
% Use [A, B, B1, D, a, b, b1] = ABfilt(...) to return
% filter coefficient vectors a, b, b1.

% Ivan Selesnick,  NYU-Poly, 2012

b1 = 1;
for i = 1:2*deg-1
    b1 = conv(b1, [-1 1]);
end
b1 = b1 * (-1)^deg;
b = conv(b1, [-1 1]);

omc = 2*pi*fc;
t = ((1-cos(omc))/(1+cos(omc)))^deg;

a = 1;
for i = 1:deg
    a = conv(a, [1 2 1]);
end
a = b + t*a;

A = spdiags( a(ones(N-2*deg,1), :), -deg:deg, N-2*deg, N-2*deg);    % A: Symmetric banded matrix
B1 = spdiags(b1(ones(N,1), :), 0:2*deg-1, N-2*deg, N-1);            % B1: banded matrix
e = ones(N, 1);
D = spdiags([-e e] , 0:1, N-1, N);
B = B1 * D;                                                         % B: banded matrix

% verify that B = B1*D
x = rand(N,1);
err = B1*diff(x) - B*x;
if max(abs(err)) > 1e-10
    disp('Error in ABfilt (B1*D not equal to B)')
end
