% Clear workspace
clear
close(findall(0,'Type','figure'));
clc;


f = @(x) (-360).*x(2)+0.45E3.*((-30)+((30+(-1).*0+x(1)).^2+((-1).*0+ ...
  x(2)).^2).^0.5E0).^2+0.3E3.*((-30)+((30+(-1).*x(1)+0).^2+((-1) ...
  .*x(2)+0).^2).^0.5E0).^2;
f = @(x) x(1) + x(1).^2+(-1).*x(2)+(-3).*x(1).*x(2)+4.*x(2).^2;
f = @(x) (11+(-1).*x(1)+(-1).*x(2)).^2+(1+x(1)+10.*x(2)+(-1).*x(1).*x(2)).^2;

f = @(x) (1 - x(1)).^2 + 100 .* (x(2) - x(1).^2).^2;
f = @(x) (x(1)).^4 + (x(2)).^4;

%f = @(x) (x(1).^2 + x(2).^2);
p = [3;2]; % [10;2]; %[0.01;-0.10] [-1;-3];

lin_name = 'Bisseção'; %Bisseção Seção Áurea
dir_name = 'MD'; %Univariante Powell SD FR NR BFGS MD
[~,p] = FindMinLinSearch(f,p,lin_name,dir_name);

% method_name = 'dogleg'; %dogleg Levenberg
% [~,p] = FindMinTrustRegion(f,p,method_name);
% p
