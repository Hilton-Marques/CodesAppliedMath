M = 3000;
data = xlsread('A��es_Tecnologia','A��es_Tecnologia');
v = data(1:10,2); %N�s e elementos

x = 0;
for i=1:numel(v)
C = nchoosek(v,i);
xc = sum(prod(C,2));
x = x + xc;
end
F = M*(1+x)
(F-M)*100/M
%Conclus�o:
%Favore�a as aplica��es que tem menores 
%oscila��es ao longo do m�s

