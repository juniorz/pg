function y = ackley(x)
% 
% Ackley function.
% Matlab Code by A. Hedar (Sep. 29, 2005).
% The number of variables n should be adjusted below.
% The default value of n=2.
% 
n=length(x);

%  s1=0; s2=0;
%  for i=1:n
%     s1 = s1 + x(i)*x(i);
%     s2 = s2 + cos(2*pi*x(i));
%  end
y = -20*exp(-0.2*sqrt((1/n)*sum(x.*x))) - exp((1/n)*sum(cos(2*pi*x))) + 20 + exp(1);



