function [out]=schaffer(a)
% f6.m
% Schaffer's F6 function
% commonly used to test optimization/global minimization problems
%
% z = 0.5+ (sin^2(sqrt(x^2+y^2))-0.5)/((1+0.01*(x^2+y^2))^2)

x = a(1);
y = a(2);


num=sin(sqrt(x^2+y^2))^2 - 0.5;
den=(1.0+0.001*(x^2+y^2))^2;

out = 0.5 + num/den;


