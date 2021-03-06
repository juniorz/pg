function [out]=rastrigin(x)
% Rastrigin.m
% Rastrigin function
% commonly used to test optimization/global minimization problems
%  out= 20 + x^2 +y^2 -10*cos(2*pi*x) -10*cos(2*pi*y);

% function y = rast(x)
% 
% Rastrigin function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 2.
% 
n = length(x); 
s = 0;
for j = 1:n
    s = s + (x(j)^2 - 10*cos(2*pi*x(j)));
end
out = 10*n + s;
