clc;
clear;
close all;

points = xlsread('Input','Planilha1');
j = 1;
Ci = cell(4,1);
for i = 1:4
    nCi = points(j,1);
    j = j+1;
    Ci{i} = points(j:j+nCi-1,:);
    j = j + nCi;
end
figure
hold on
for i = 1:4
    c = Ci{i};
    plot(c(:,1),c(:,2),'o-');
end
%% Correção do input
C1 = Ci{1};
C2 = Ci{2};
C3 = Ci{3};
C4 = Ci{4};
nC1 = size(C1,1);
nC2 = size(C2,1);
nC3 = size(C3,1);
nC4 = size(C4,1);
if nC1 > nC2
    Ci{2} = correction(C2,C1);
elseif nC2 > nC1
    Ci{1} = correction(C1,C2);
end
if nC3 > nC4
    Ci{4} = correction(C4,C3);
elseif nC4 > nC3
    Ci{3} = correction(C3,C4);
end
figure
hold on
for i = 1:4
    c = Ci{i};
    plot(c(:,1),c(:,2),'o-');
end
%% Start
nu = size(Ci{1},1);
nv = size(Ci{3},1);
u = linspace(0,1,nu);
v = linspace(0,1,nv);
pvx = (1-v)'.*Ci{1}(:,1)' + v'.*Ci{2}(:,1)';
pvy = (1-v)'.*Ci{1}(:,2)' + v'.*Ci{2}(:,2)';
pux = (1-u).*Ci{3}(:,1) + u.*Ci{4}(:,1);
puy = (1-u).*Ci{3}(:,2) + u.*Ci{4}(:,2);
ppx = (1 - u) .* (1 - v)' * Ci{1}(1,1) + ...
    u.*(1 - v)'*Ci{1}(end,1) + u.*v'*Ci{2}(end,1) + ...
    (1-u).* v' * Ci{2}(1,1);
ppy = (1 - u) .* (1 - v)' * Ci{1}(1,2) + ...
    u.*(1 - v)'*Ci{1}(end,2) + u.*v'*Ci{2}(end,2) + ...
    (1-u).* v' * Ci{2}(1,2);

meshx = pvx + pux  - ppx ;
meshy = pvy + puy - ppy ;
fileID = fopen('u.txt','w');
fprintf(fileID,'%10s %10s %10s %10s %10s %10s %10s %10s\n',...
    'x1','x2','x3','x4','y1','y2','y3','y4');
for i = 1:nv-1
    for j = 1: nu-1
        p1 = [meshx(i,j),meshy(i,j)];
        p2 = [meshx(i,j+1),meshy(i,j+1)];
        p3 = [meshx(i+1,j+1),meshy(i+1,j+1)];
        p4 = [meshx(i+1,j),meshy(i+1,j)];
        px = [p1(1),p2(1),p3(1),p4(1)];
        py = [p1(2),p2(2),p3(2),p4(2)];
        fprintf(fileID,'%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n',...
            [px,py]);
        drawRect(p1,p2,p3,p4);
    end
end
 fclose(fileID);


function newCMenor = correction(c_menor,c_maior)
nc_menor = size(c_menor,1);
nc_maior = size(c_maior,1);
newCMenor = zeros(nc_maior,2);
newCMenor(1,:) = c_menor(1,:);
newCMenor(end,:) = c_menor(end,:);
u = 0:1/(nc_maior-1):1;
h = 1/(nc_menor-1);
for i = 2 : nc_maior-1
    t = u(i);
    pos = fix(t/h);
    a = c_menor(pos + 1,:);
    b = c_menor(pos + 2,:);
    p = a + ((t-pos*h)/h)*(b-a);
    newCMenor(i,:) = p;
end
end
function drawRect(p1,p2,p3,p4)
  line([p1(1),p2(1)],[p1(2),p2(2)]);
  line([p2(1),p3(1)],[p2(2),p3(2)]);
  line([p3(1),p4(1)],[p3(2),p4(2)]);
  line([p4(1),p1(1)],[p4(2),p1(2)]);
end