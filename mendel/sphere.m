function [out]=sphere(x)
% 
% Sphere function 
% Matlab Code by A. Hedar (Nov. 23, 2005).
% The number of variables n should be adjusted below.
% The default value of n = 30.
% 

%  s = 0;
%  for j = 1:length(x)
%      s = s + x(j)^2; 
%  end
%  out = s;

out = sum(x.*x);
