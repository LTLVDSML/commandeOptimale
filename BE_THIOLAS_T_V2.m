function[] = BE_THIOLAS_T_V2()

clc 
clear all
close all

A = [-0.746 0.006 -1 0.0369;     %matrice de commande
    -12.9 -0.746 0.387 0;
    4.31 0.024 -0.174 0;
    0 1 0 0];

B = [0.0012 0.0092;            %matrice observation
    6.05 0.952;
    -0.416 -1.76;
    0 0];

C = [20 0 0 0;           %matrice regroupant les coeffs 
    0 0 0 0;             %de la derivation de la fonction par 
    0 0 0 0;             %les elements de X
    0 0 0 20];


%%%%%%%%%%%%%%%%%%%%Méthode de controle optimal%%%%%%%%%%%%%%%%%%%%%%%%%
global G

G = [A (-B*B');  %on crée la matrice G pour ode
     -C -A'];

T = 4;           %fin de la mesure
dt = 0.1;        %pas de temps
t = 0:dt:T;

global xo

%conditions a t=0
xo = [0.087266;0;0;0.122173];  %x a t = 0
p = [0;0;0;0];                %p a t = 0


sol = bvp4c(@ode,@bc,bvpinit(t,zeros(8,1)));  %on determine les valeurs optimale dans le temps grace a bvp4c   

figure(1)                  %on affiche les commandes optimales pour chaque parametre
subplot(2,2,1)
plot(sol.x,sol.y(1,:))
title('beta(t)')
subplot(2,2,2)
plot(sol.x,sol.y(2,:))
title('p(t)')
subplot(2,2,3)
plot(sol.x,sol.y(3,:))
title('r(t)')
subplot(2,2,4)
plot(sol.x,sol.y(4,:))
title('phi(t)')

figure(2)                      %on affiche les couts de chaque parametre
subplot(2,2,1)
plot(sol.x,sol.y(5,:))
title('p(beta(t))')
subplot(2,2,2)
plot(sol.x,sol.y(6,:))
title('p(p(t))')
subplot(2,2,3)
plot(sol.x,sol.y(7,:))
title('p(r(t))')
subplot(2,2,4)
plot(sol.x,sol.y(8,:))
title('p(phi(t))')

%%%%%%%%%%%%%%%%%%%%%%%%%%Méthode directe%%%%%%%%%%%%%%%%%%%%%%

T = 4;           %fin de la mesure
dt = 0.01;        %pas de temps (ne prendre que des diviseurs de 4)
t = 0:dt:T;

Q = zeros(6*(T/dt)+4,6*(T/dt)+6);  %matrice Q initialisation
Qk = zeros(6,6);                   %initialisation matrice Qk
Qk(1,1) = 20*dt;                   %qui sera utilisee pour former 
Qk(4,4) = 20*dt;                   %la matrice Q
Qk(5,5) = dt;
Qk(6,6) = dt;

Q(1:6,1:6) = Qk;        %creation matrice Q

for i = 1:1+T/dt-1      %on copie la matrice Qk toutes les 6 lignes et colonnes
                        %pour former la matrice Q
    Q(1+i*6,1+i*6) = Qk(1,1);
    Q(4+i*6,4+i*6) = Qk(4,4);
    Q(5+i*6,5+i*6) = Qk(5,5);
    Q(6+i*6,6+i*6) = Qk(6,6);
    
end


H = zeros(4*(T/dt+1),6*(T/dt+1));   %matrice H initialisation

[ligH,colH] = size(H);       %generalisation
H(1:4,1:4) = eye(4);         %matrice identite 4 premieres lignes et colonnes de H
H(ligH-3:ligH,colH-5:colH-2) = eye(4);  %matrice identite 4 dernieres ligne de H

Hk = zeros(4,10);               %initialisation matrice Hk
Hk(1:4,1:4) = eye(4)+dt.*A;     %qui sera utilisee pour former
Hk(1:4,5:6) = dt.*B;            %la matrice H
Hk(1:4,7:10) = -eye(4);

for i = 1:T/dt-1    %on copie la matrice Hk toutes les 6 lignes et colonnes
                    %pour former la matrice H
    H(1+i*4:4+i*4,1+(i-1)*6:(i-1)*6+10) = Hk;
 
end 

g = zeros(4*(T/dt+1),1);  %initialisation matrice contrainte
g(1:4) = xo;              %conditions initiales

X = H\g;                 %matrice comprenant les valeurs des variables d etat au cours du temmps
a = length(X);

p = zeros(4*(T/dt+1),1);   %initialisation matrice couts au cours du temps
b = length(p);

x = [X;p];                %matrice regroupant valeurs de X et couts au cours du temps

%methode de Newton

tol = 0.0001;                %tolerance

while norm(gradient(Q,X,H,g,p),2)>tol     

x = x-hessienne(Q,H)\gradient(Q,X,H,g,p);
X = x(1:a);
p = x(a+1:a+b);
end

x1 = zeros(1,length(X)/6);   %on intialise les vecteurs regroupant les valeurs de 
x2 = zeros(1,length(X)/6);   %chaque variable d etat
x3 = zeros(1,length(X)/6);
x4 = zeros(1,length(X)/6);

p1 = zeros(1,length(X)/6);   %on intialise les vecteurs regroupant les valeurs de
p2 = zeros(1,length(X)/6);   %chaque variable d etat
p3 = zeros(1,length(X)/6);
p4 = zeros(1,length(X)/6);

for i=1:length(X)/6
    x1(i) = X(1+6*(i-1));   %on separe les composantes de la matrice X
    x2(i) = X(2+6*(i-1));   %pour obtenir les valeurs des varaibles d etat 
    x3(i) = X(3+6*(i-1));   %dans le temps
    x4(i) = X(4+6*(i-1));

    p1(i) = p(1+4*(i-1));   %on separe les composantes de la matrice p
    p2(i) = p(2+4*(i-1));   %pour obtenir les valeurs des couts
    p3(i) = p(3+4*(i-1));   %dans le temps
    p4(i) = p(4+4*(i-1));
end

figure(3)         %on affiche les variables d etat
subplot(2,2,1)
plot(t,x1)
title('beta(t)')
subplot(2,2,2)
plot(t,x2)
title('p(t)')
subplot(2,2,3)
plot(t,x3)
title('r(t)')
subplot(2,2,4)
plot(t,x4)
title('phi(t)')

figure(4)        %on affiche les couts associes
subplot(2,2,1)
plot(t,p1)
title('p(beta(t))')
subplot(2,2,2)
plot(t,p2)
title('p(p(t))')
subplot(2,2,3)
plot(t,p3)
title('p(r(t))')
subplot(2,2,4)
plot(t,p4)
title('p(phi(t))')

end


%REMARQUE :
%On observe que la fonction bvp4c permet d otenir un resultat avec une 
%precision similaire a la methode directe tout en necessitant un pas de 
%temps moins petit
%Cette fonction permet donc d obtenir des resultats aussi precis que la methode directe
%tout en necessitant moins de ressources


%%%%%%%%%%%%%%%%%%%%%%%%%Fonctions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù%%

%%%%%%%Methode optimale%%%%%%%%%%%%%%

function od = ode(x,y)
global G
od = G*y;   %on caclule la derivee de y

end

function psi = bc(y0,yT)  %nul si atteint les limites

global xo
x0 = y0(1:4);    %on retire les elements de p
xT = yT(1:4);

  psi = [x0-xo ;  
      xT];                          
              %on verifie les conditions aux limites
end 


%%%%%%%Methode directe%%%%%%%%%%%%%

function gr = gradient(Q,X,H,g,p) %matrice gradient 

gr = [Q*X+H'*p;   
    H*X-g];

end

function he = hessienne(Q,H)   %matrice hessienne

[ligH,colH] = size(H);      %generalisation

he = [Q H';
      H zeros(ligH,ligH)];

end
