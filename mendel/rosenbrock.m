function [out]=rosenbrock(x)

n = length(x);

sum = 0;
for i = 1:n-1;
   sum = sum + 100*(x(i+1)-x(i)^2)^2 + (x(i)-1)^2;
end
out = sum;
