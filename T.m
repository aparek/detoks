function y = T(x)
% function y = T(x)
% 
% Applies Teager Kaiser Operator to input x
% y - output
%
% Ankit Parekh
% LAST EDIT: 7/21/2018
% ankit.parekh@mssm.edu

for i = 2:length(x)-1
    y(i) = x(i)^2 - x(i-1)*x(i+1);
end

end

