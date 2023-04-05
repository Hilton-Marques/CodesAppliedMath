close all
clc;
clear;

F = eye(7);
PP = [[10011.         0.         0.         0.     10000.         0.      0.    ];...
 [    0.     10011.         0.         0.         0.     10000.      0.    ];...
 [    0.         0.     10011.         0.         0.         0.  10000.    ];...
 [    0.         0.         0.        11.         0.         0.      0.    ];...
 [10000.         0.         0.         0.     10000.01       0.      0.    ];...
 [    0.     10000.         0.         0.         0.     10000.01      0.    ];...
 [    0.         0.     10000.         0.         0.         0.  10000.0001]];
F(1,5) = 1;
F(2,6) = 1;
F(3,7) = 1; % mudança de estado;
H = eye(4,7); % medições pela detecção 

% Recall : z = H*x + R , as medições são u,v,s,r posição,área e razão de
% aspecto
R = eye(4,4); 
R(3:4,3:4) = 10*R(3:4,3:4); % maior incerteza para a área e para a razão de aspecto
R
% Incerteza inicial P
P = eye(7);
P = 10*P;
P(5:end,5:end) = 1000*P(5:end,5:end); %penalize as variaveis nao observaveis (todas as velocidades)
P
% a nova variancia pela predição sera PFP' + Q
Q = 10*eye(7);
Q(5:end,5:end) = 0.01*Q(5:end,5:end); %penalize as variaveis nao observaveis (todas as velocidades)
Q(end,end) = 0.01*Q(end,end);
Q
%teste 
x  = [3.84660000e+02; 2.01624900e+02; 2.14114435e+04; 3.58902764e-01; 0.00000000e+00; 0.00000000e+00; 0.00000000e+00];
x = [ 156;    282;    61008;  0.44086;      0;     0;      0];
xpred = F*x 
Ppred = F*P*F' + Q
z = [3.92976000e+02 ; 2.09915000e+02; 2.08186462e+04; 3.57012117e-01];
z = [154;281; 59200; 0.43243];
i = (z - H*xpred)
K = Ppred*H'*inv(H*Ppred*H' + R);
v = H'*inv(H*Ppred*H' + R)*i;
Ppred(6,:)

xup = xpred + K*i
Pup = Ppred - K*H*Ppred
