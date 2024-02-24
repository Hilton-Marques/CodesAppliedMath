% Uma massa simples
t = 0.1;
A = [1 0; t/2 1];
B = [(1 - (t^2)/2) t;-t/2 1];
M = (inv(A)*B);

% Duas massas
t = .64/1;%Teste que com t = 0.63 o resultado está proximo da instabilidade
          %Com t = 0.64 o resultado explode
M=[9 0; 0 1];  K = [81 -6; -6 6] ;
uold = [1;0]; vold = [0;0];
for i = 1:100
    vh = vold - (t/2)*(inv(M)*(K*uold));
    unew = uold + t*vh;
    vnew = vh - (t/2)*(inv(M)*(K*unew));
    uold = unew;
    vold = vnew;
end
% u de t = 63
unew
