M = 3000;
data = xlsread('Ações_Tecnologia','Ações_Tecnologia');
v = data(1:10,2); %Nós e elementos

x = 0;
for i=1:numel(v)
C = nchoosek(v,i);
xc = sum(prod(C,2));
x = x + xc;
end
F = M*(1+x)
(F-M)*100/M
%Conclusão:
%Favoreça as aplicações que tem menores 
%oscilações ao longo do mês

