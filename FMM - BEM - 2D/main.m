% Clear workspace
clear
close(findall(0,'Type','figure'));
%clc
%Input data
points = xlsread('Input','Planilha1');
nEl = points(1,2);
nFaces = nEl;
pointsCor = points(2:end,1:2)';
typeBC = points(2:end,4);
BC = points(2:end,5);
%q = zeros(2,nFaces);
nNos = 1;
nodeInt = 1;
nNode = nEl;
intCord = zeros(2,nodeInt);
Elements(nEl) = Element();
Nodes(nNode) = Node();
for i = 1:nEl
    q(1,1) = typeBC(i);
    q(2,1) = BC(i);
    face = Face(q);
    nelCord(1:2,1) = pointsCor(:,i);
    if (i == nEl)
        nelCord(3:4,1) = pointsCor(:,1);
        r = (nelCord(3:4,1) - nelCord(1:2,1));
        len = (dot(r,r))^0.5;
        Elements(i) = Element(nNos,len,face,nelCord,i);
        continue;
    end
    nelCord(3:4,1) = pointsCor(:,i+1);
    r = (nelCord(3:4,1) - nelCord(1:2,1));
    len = (dot(r,r))^0.5;
    Elements(i) = Element(nNos,len,face,nelCord,i);
end
for i = 1:nNode
    elPosition = Elements(i).position;
    pos = (elPosition(1:2) + elPosition(3:4))/2;
    Nodes(i) = Node(pos,Elements(i));
end
hold on
DrawMesh(Nodes,Elements)
tic
k = 15; %truncation of FMM
capacity = 5; % max capacity for each leaf
Field(Nodes,capacity,k);
toc

function DrawMesh(Nodes,Elements)
% Dúvida : como acessar todos os elementos de um objeto sem
% precisar do for??
fac1 = 1 + length(Nodes)/40;
fac2 = 1 + (1-1)/10;

for i = 1 : length(Nodes)
    xn = Nodes(i).position(1,1);
    yn = Nodes(i).position(2,1);
    circle(xn,yn,(1/40)/(fac1*fac2),[0 0 1]);
    text(xn,yn,num2str(i),'HorizontalAlignment','right')
end
for i = 1 : length(Elements)
    xn = Elements(i).position(1,1);
    yn = Elements(i).position(2,1);
    xn2 = Elements(i).position(3,1);
    yn2 = Elements(i).position(4,1);
    line([xn-((0.03)/fac1), xn+(0.03)/fac1], [(yn-(0.03)/fac1),(yn +(0.03)/fac1)], 'Color', [0 0 0]);
    line([xn,xn2],[yn,yn2],'Color', [0 0 0]);
    %line([this.els(1).cord(1,1),this.els(end).cord(1,1)],[this.els(1).cord(2,1),this.els(end).cord(2,1)],'Color', [0 0 0]);
end
end
function circle(x,y,r,c)
circ = 0 : pi/50 : 2*pi;
xcirc = x + r * cos(circ);
ycirc = y + r * sin(circ);
plot(xcirc, ycirc, 'color', c);
end





