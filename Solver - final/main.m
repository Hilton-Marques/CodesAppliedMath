%% Clear memory
clc;
clear; 
close all;

%% Test model
nfaces = 4;                    % number of faces from geometry
nelperface = 1;                % number of elements per face
nnel = nelperface*nfaces;       %numer of elements
no = 3;                        % nodes per element 
continum = true;               % flag for continuum nodes

if continum == true
    nnode = (no-1)*nnel;              %number of nodes
else
    nnode = (no)*nnel;
end
               %number of int nodes

facescord = zeros(nfaces,4);    % cordinates of faces
cond = ones(1,nfaces);             % Conductivity per face
q = zeros(2,nfaces);             % array of prescibed contourn values [1] para gradiente conhecido [0] for temperatura conhecida
nodescor = zeros(2,nnode);      % array for nodes coordinates [x] and [y]
nelface = zeros(1,nnel);        % face and boundary condition for each element
nelcord = zeros(4,nnel);        % Coord of begin and start for each element
nodeelement = zeros(2,nnode);   % which element node belong
elementnode = zeros(no,nnel);
cond =1;
facescord(1,:) = [0 0 1 0];
facescord(2,:) = [1 0 1 1];
facescord(3,:) = [1 1 0 1];
facescord(4,:) = [0 1 0 0];

q(1,1) = 1;                 % Boundary condition per face ...                          
q(2,1) = 0;
q(1,2) = 0;
q(2,2) = 200;
q(1,3) = 1;
q(2,3) = 0;
q(1,4) = 0;
q(2,4) = 100;

% obtain the coordinates per elements
k=1;
for i =1:nfaces
   
    vectorface = [(facescord(i,3)-facescord(i,1)) ;  (facescord(i,4)-facescord(i,2))];
    teta1 = acos(vectorface(1)/norm(vectorface));
    teta2 = asin(vectorface(2)/norm(vectorface));
    for j=1:nelperface
        nelcord(3,k) = nelcord(1,k) + cos(teta1)*(1/nelperface)*norm(vectorface);
        nelcord(4,k) = nelcord(2,k) + sin(teta2)*(1/nelperface)*norm(vectorface);
        nelface(1,k) = i;
        if k > nelperface*nfaces - 1
            break
        end
        nelcord(1,k+1) = nelcord(3,k); 
        nelcord(2,k+1) = nelcord(4,k);        
        k = k+1;
         
    end
end

c =1;
if continum == false
    for i = 1: nnel
    vectorelement = [(nelcord(3,i)-nelcord(1,i)) ;  (nelcord(4,i)-nelcord(2,i))];
    teta1 = acos(vectorelement(1)/norm(vectorelement));
    teta2 = asin(vectorelement(2)/norm(vectorelement));
        for j=1:no
        nodescor(1,c) = nelcord(1,i) + (j)*cos(teta1)*(1/(no+1))*norm(vectorelement);
        nodescor(2,c) = nelcord(2,i) + (j)*sin(teta2)*(1/(no+1))*norm(vectorelement); 
        nodeelement(1,c) = i;
        nodeelement(2,c) = i;
        c = c+1;
        end     

    end
else
    for i = 1: nnel
        
        vectorelement = [(nelcord(3,i)-nelcord(1,i)) ;  (nelcord(4,i)-nelcord(2,i))];
        teta1 = acos(vectorelement(1)/norm(vectorelement));
        teta2 = asin(vectorelement(2)/norm(vectorelement));
        nodescor(1,c) = nelcord(1,i);
        nodescor(2,c) = nelcord(2,i);
        if c > 1
            nodeelement(2,c) = i;
        end
        for j=1:no-1
            if c == nnode
                nodeelement(1,1) = i;
                nodeelement(2,nnode) = i;
                elementnode(no,i) = 1;
                elementnode(1,i) = no*(nnel -1) - (nnel -2);
                break
            end
            nodescor(1,c+1) = nelcord(1,i) + (j)*cos(teta1)*(1/(no-1))*norm(vectorelement);
            nodescor(2,c+1) = nelcord(2,i) + (j)*sin(teta2)*(1/(no-1))*norm(vectorelement);
            nodeelement(2,c) = i;
            nodeelement(1,c+1) = i;
            elementnode(j,i) = c;
            elementnode(j+1,i) = c+1;
            c = c + 1;
        end
        
    end
end
nodeelement([2 1],:) =  nodeelement; %deixando na ordem das funções de forma e.g F2*el2 + F1*el1

x = linspace(0.001,0.999,10);
y = linspace(0.99999,0.99999,10);
ym = linspace(0.00001,0.00001,10);
nodeint = 2*numel(x);    
intcord = zeros(2,nodeint);     % coordinates for each int nod
intcord = [ x x; y ym];

%% Solve model iteratively and print each stage
potential = PotentialSolver(nfaces,no,nnel,nodescor,cond,nelface,nelcord, ...
    nodeint, intcord,facescord,q,nodeelement,nnode,elementnode);

potential.Solve();
potential.Solveint();
potential.DrawMesh();

