clc        %%%%%%%Page 2.3 du cours%%%%%%%%%
clear all 
close all

x = [0 0 0]';
tol = 0.00000001;
i = 0;

while norm(gradient(x),2)>tol

x = x-hessienne(x)\gradient(x);
i = i+1;
end

X = x(1)
U = x(2)
lambda = x(3)

iterations =i
%fonctions


function dg = gradient(x)
x1=x(1);
x2=x(2);
x3=x(3);
a=1;
b=2;
c=0.5;
m=4;
dg = zeros(3,1);
dg(1) = x1/a^2+x3;
dg(2) = x2/b^2+m*x3;
dg(3) = x1+m*x2-c;
end

function dg2 = hessienne(x)
x1=x(1);
x2=x(2);
x3=x(3);
a=1;
b=2;
c=0.5;
m=4;
dg2 = zeros(3,3);
dg2(1,1) = 1/a^2;
dg2(1,2) = 0;
dg2(1,3) = 1; 
dg2(2,1) = 0;
dg2(2,2) = 1/b^2;
dg2(2,3) = m;
dg2(3,1) = 1;
dg2(3,2) = m;
dg2(3,3) = 0;
end