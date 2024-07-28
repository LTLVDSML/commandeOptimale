clc        %%%%%%%Page 3.5 du cours%%%%%%%%%
clear all 
close all

x = 2;
lambda = [1 1]';
tol = 0.00000001;
rho = 0.1;
i = 0;
ii = 0;
%C(x) = [x-1; -x]

while norm(lambda.*[x-1; -x],2)>tol || max([x-1; -x])>tol

    while norm(gradient(x,lambda),2)>tol
        i =0;
        x = x -hessienne(x,lambda)\gradient(x,lambda);
        i = i+1;
        
    end
    
lambda = lambda+rho*[x-1; -x];
lambda(lambda<0) = 0;
ii = ii+1;
end

x
lambda

%minimum en x = 0, lambda = [0,1]

i

ii


%fonctions


function dg = gradient(x,lambda)
x1=x;
lambda1 = lambda(1);
lambda2 = lambda(2);
dg = x1+1+lambda1-lambda2;
end

function dg2 = hessienne(x,lambda)
x1=x;
dg2 = 1;
end