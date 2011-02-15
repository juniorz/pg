function[z]=griewank(x)
%x1 = x(1);
%x2 = x(2);


cosprod = 1;
quadsum = sum(x .* x);
for i=1:length(x)
    cosprod = cosprod * cos(x(i)/sqrt(i));
end
z = quadsum/4000 - cosprod + 1;

% Griewank.m -Geenralized Rastrigin
%-600<=x(i)<= 600 
%y = cos(x(1)/sqrt(1))*cos(x(2)/sqrt(2))
%a = x(2)/sqrt(2)

%z = (x1^2 + x2^2)/4000 - cos(x1/sqrt(1))*cos(x2/sqrt(2)) + 1;




