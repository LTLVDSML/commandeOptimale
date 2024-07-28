clc         %%%%%%%Page 1.4 du cours%%%%%%%%%
clear all
close all

%Methode grad a pas csts

x = [0 0]';
rho = 0.001;
tol = 0.0000001;

while norm(gradient(x),2)>tol
    
    x = x-rho*gradient(x);
    
end

x

%methode de Newton

x = [0 0]';

while norm(gradient(x),2)>tol

x = x-hessienne(x)\gradient(x);

end

x

%minimum en [1 1]


%fonctions

function dg = gradient(x)
x1=x(1);
x2=x(2);
dg = zeros(2,1);
dg(1) = 400*x1^3 +2*x1 - 400*x1*x2 - 2;
dg(2) = 200*x2 - 200*x1^2;
end

function dg2 = hessienne(x)
x1=x(1);
x2=x(2);
dg2 = zeros(2,2);
dg2(1,1) = -400*x2 +2 +1200*x1^2;
dg2(1,2) = -400*x1;
dg2(2,1) = -400*x1;
dg2(2,2) = 200;
end