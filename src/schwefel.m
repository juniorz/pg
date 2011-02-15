function y = schwefel(x)
% 
% Schwefel function
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 2.
% 
n = length(x);
i=1:n;
%  s = sum(-(x+1.105*i).*sin(sqrt(abs((x+1.105*i)))));
s = -sum(x.*sin(sqrt(abs(x))));
%y = s;
y = 418.9829*n + s;
