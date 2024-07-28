%function[] = BE_THIOLAS_T()

clc 
clear all
close all

A = [-0.746 0.006 -1 0.0369;
    -12.9 -0.746 0.387 0;
    4.31 0.024 -0.174 0;
    0 1 0 0];

B = [0.0012 0.0092;
    6.05 0.952;
    -0.416 -1.76;
    0 0];

C = [20 0 0 0;
    0 0 0 0;
    0 0 0 0;
    0 0 0 20];


%%%%%%%%%%%%%%%%%%%%Méthode de controle optimal%%%%%%%%%%%%%%%%%%%%%%%%%
global G

G = [A (-B*B');  %on globalise la matrice G pour ode
     -C -A'];

T = 4;
dt = 0.1;
t = 0:dt:T;

global xo

xo = [0.087266;0;0;0.122173];  %x a t = 0
p = [0;0;0;0];                %p a t = 0
% y = [x; p];


sol = bvp4c(@ode,@bc,bvpinit(t,zeros(8,1)));

figure(1)
subplot(2,2,1)
plot(sol.x,sol.y(1,:))
title('x1')
subplot(2,2,2)
plot(sol.x,sol.y(2,:))
title('x2')
subplot(2,2,3)
plot(sol.x,sol.y(3,:))
title('x3')
subplot(2,2,4)
plot(sol.x,sol.y(4,:))
title('x4')

figure(2)
subplot(2,2,1)
plot(sol.x,sol.y(5,:))
title('p1')
subplot(2,2,2)
plot(sol.x,sol.y(6,:))
title('p2')
subplot(2,2,3)
plot(sol.x,sol.y(7,:))
title('p3')
subplot(2,2,4)
plot(sol.x,sol.y(8,:))
title('p4')

%%%%%%%%%%%%%%%%%%%%%%%%%%Méthode directe%%%%%%%%%%%%%%%%%%%%%

T = 4;
dt = 0.1;

Q = zeros(6*T/dt,6*T/dt);
Qk = zeros(6,6);
Qk(1,1) = 10;
Qk(4,4) = 10;
Qk(5,5) = 0.5;
Qk(6,6) = 0.5;

Q(1:6,1:6) = Qk;

for i = 1:1+T/dt-2
    
    Q(1+i*6,1+i*6) = Qk(1,1);
    Q(4+i*6,4+i*6) = Qk(4,4);
    Q(5+i*6,5+i*6) = Qk(5,5);
    Q(6+i*6,6+i*6) = Qk(6,6);
    
end

H = zeros(6*T/dt,6*T/dt);
H(1:6,1:6) = eye(6);

Hk = zeros(6,6);
Hk(1:4,1:4) = eye(4)+A;
Hk(1:4,5:6) = B;

for i = 1:T/dt-2
    
    H(1+i*6:6+i*6,1+i*6:6+i*6) = Hk;
 
end

uo = [0;0];
x = xo;
X = [x;uo];

g = zeros(6*T/dt,1);
g(1:6) = X;

for i = 1:T/dt-2
    
    X = X + dt*Hk*X;
    g(i*6+1:i*6+6) = X;
    
end

p = zeros(1,6*T/dt);

%methode de Newton
tol = 0.00001;
X = zeros(240,1);

while norm(gradient(Q,X,H,g,p),2)>tol

x = x-hessienne(Q,H)\gradient(Q,X,H,g,p);

end

%end


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

function gr = gradient(Q,X,H,g,p)

gr = [Q*X+H*p';
    H*X-g];

end

function he = hessienne(Q,H)

he = [Q H;
      H zeros(240,240)];

end
