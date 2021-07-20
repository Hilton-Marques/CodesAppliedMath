% Clear workspace
clear
close(findall(0,'Type','figure'));
clc;
% Choice draw or input
flag = 1;
if flag
data = xlsread('Input','Planilha1');
nPoints = data(1,2);
nEl = data(1,3);
nFaces = nEl;
nNodes = nEl;
nNos = 1;
Faces(nFaces) = Face();
Elements(nEl) = Element();
Nodes(nEl) = Node();
points = data(2:nPoints+1,1:3);
startEl = nPoints+1;
for i = 1:nFaces
    Faces(i).q = data(startEl+i,4:5);
end
for i = 1:nEl
    Elements(i) = Element(Faces(i),nNos,data(startEl+i,1:3),points,i);
end
intNode(1) = Node;
else
points = [1,0,0;...
    1,0,1;...
    1,1,1;...
    1,1,0;...
    0,0,0;...
    0,0,1;...
    0,1,1;...
    0,1,0];
nFaces = 6;
nEl = 12;
nNodes = 12;
nNos = 1; % Number of nodes per element
%% Boundary Conditions
Faces(nFaces) = Face();
Faces(1).q = [1;1];
Faces(2).q = [1;2];
Faces(3).q = [1;3];
Faces(4).q = [1;4];
Faces(5).q = [1;5];
Faces(6).q = [1;6];
%% Elements
Elements(nEl) = Element();
Elements(1) = Element(Faces(1),nNos,[1,4,3],points,1);
Elements(2) = Element(Faces(1),nNos,[1,3,2],points,2);
Elements(3) = Element(Faces(2),nNos,[4,8,7],points,3);
Elements(4) = Element(Faces(2),nNos,[4,7,3],points,4);
Elements(5) = Element(Faces(3),nNos,[8,5,6],points,5);
Elements(6) = Element(Faces(3),nNos,[8,6,7],points,6);
Elements(7) = Element(Faces(4),nNos,[5,1,2],points,7);
Elements(8) = Element(Faces(4),nNos,[5,2,6],points,8);
Elements(9) = Element(Faces(5),nNos,[4,5,8],points,9);
Elements(10) = Element(Faces(5),nNos,[4,1,5],points,10);
Elements(11) = Element(Faces(6),nNos,[6,2,3],points,11);
Elements(12) = Element(Faces(6),nNos,[6,3,7],points,12);
end

%% Nodes

Nodes(nNodes) = Node();
for i = 1:nNodes
    Nodes(i) = Node(Elements(i));
end

%% FMM
n = 3; %truncation of FMM
capacity = 40; % max capacity for each leaf
ng = 3; %Gauss points for integration
tic
fmm = Field(Nodes,capacity,n,ng);
fmm.solver.Ax
toc
keyboard
%% Draw Mesh
hold on
Draw(Nodes,nNodes,Elements,nEl)
function Draw(Nodes,nNodes,Elements,nEl)
view(-34,31);
for i = 1:nNodes
    p = Nodes(i).pos;
    plot3(p(1),p(2),p(3),'o','Color','Red');
    text(p(1),p(2),p(3),num2str(i),'HorizontalAlignment','left');
end
for i = 1:nEl
    element = Elements(i);
    element.plot;
end
end

