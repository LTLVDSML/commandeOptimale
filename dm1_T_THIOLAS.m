function[]=dm1_T_THIOLAS()

clc
close all
clear all

%%%%%%%%%%%%%%%%Optimisation sans contrainte%%%%%%%%%

%Methode grad a pas csts

x = [10 10 10]';  %initialisation de x 
rho = 0.001;       %pas 
tol = 0.00000001;   %tolerance suffisamment grande pour etre le plus proche du resultat

while norm(gradient(x),2)>tol    %tant que la norme du gradient est superieure a la tolerance
    
    x = x-rho*gradient(x);    %on déplace x selon le gradient de x et rho
    
end

x

%minimum en [1 0 0]
%On donne la valeur atteinte sur le minimum
minimum = 1.5*x(1)^2+2*x(2)^2+1.5*x(3)^2+x(1)*x(3)+2*x(2)*x(3)-3*x(1)-x(3)

%methode de Newton

x = [10 10 10]'; %initialisation de x

while norm(gradient(x),2)>tol   %tant que la norme du gradient est superieure a la tolerance

x = x-hessienne(x)\gradient(x);   %on deplace x selon la division de sa hessienne par son gradient

end

x

%minimum en [1 0 0]
%On donne la valeur atteinte sur le minimum
minimum = 1.5*x(1)^2+2*x(2)^2+1.5*x(3)^2+x(1)*x(3)+2*x(2)*x(3)-3*x(1)-x(3)

%%%%%%%%%%%%%%%%Optimisation avec contraintes%%%%%%%%%

%Methode d Uzawa

x = [0.1 0.1 0.1]';   %initialisation de x
lambda = 1;          %initialisation de lambda
tol = 0.00000001;   %tolerance
rho = 0.01;       %pas
i = 0;
ii = 0;
C = [10*x(1)^2+10*x(2)^2+10*x(3)^2-1];    %contrainte


while norm(lambda*C,2)>tol || max(C)>tol
   %tant que la norme du produit de lambda par C(x) ou C(x) est superieur a la tolerance
   
    while norm(gradient2(x,lambda),2)>tol
        %tant que la norme du gradient est superieure a la tolerance
        i =0;
        x = x -hessienne2(x,lambda)\gradient2(x,lambda);  %Newton
        %x = x-rho*gradient2(x,lambda);            %grad a pas constants
        i = i+1; 
    end
    
C = [10*x(1)^2+10*x(2)^2+10*x(3)^2-1];%On redefinit la contrainte  

lambda = lambda+rho*[10*x(1)^2+10*x(2)^2+10*x(3)^2-1]; 
%on modifie lambda en lui ajoutant le produit de rho*C(x)
lambda(lambda<0) = 0;
%si lambda est négatif, on lui affecte 0
ii = ii+1;
end

x
lambda

%on verifie que la contrainte est bien respectee
condition = 10*(x(1)^2+x(2)^2+x(3)^2)-1
%On donne la valeur atteinte sur le minimum
minimum = 1.5*x(1)^2+2*x(2)^2+1.5*x(3)^2+x(1)*x(3)+2*x(2)*x(3)-3*x(1)-x(3)

%minimum en [0.3067 -0.0144 0.0757]

end

%fonctions
%%%%%%%%%%%%%%ex1%%%%%%%%%%%%%%%%
function dg = gradient(x)
x1=x(1);
x2=x(2);
x3=x(3);
dg = zeros(3,1);
dg(1) = 3*x1+x3-3;
dg(2) = 4*x2+2*x3;
dg(3) = x1+2*x2+3*x3-1;
end

function dg2 = hessienne(x)
x1=x(1);
x2=x(2);
x3=x(3);
dg2 = zeros(3,3);
dg2(1,1) = 3;
dg2(1,2) = 0;
dg2(1,3) = 1;
dg2(2,1) = 0;
dg2(2,2) = 4;
dg2(2,3) = 2;
dg2(3,1) = 1;
dg2(3,2) = 2;
dg2(3,3) = 3;
end

%%%%%%%%%%%%%%ex2%%%%%%%%%%%%%%%%
function dg = gradient2(x,lambda)
x1=x(1);
x2=x(2);
x3=x(3);
dg = zeros(3,1);
dg(1) = 3*x1+x3-3+20*lambda*x1;
dg(2) = 4*x2+2*x3+20*lambda*x2;
dg(3) = x1+2*x2+3*x3-1+20*lambda*x3;
end

function dg2 = hessienne2(x,lambda)
x1=x(1);
x2=x(2);
x3=x(3);
dg2 = zeros(3,3);
dg2(1,1) = 3+20*lambda;
dg2(1,2) = 0;
dg2(1,3) = 1;
dg2(2,1) = 0;
dg2(2,2) = 4+20*lambda;
dg2(2,3) = 2;
dg2(3,1) = 1;
dg2(3,2) = 2;
dg2(3,3) = 3+20*lambda;
end
