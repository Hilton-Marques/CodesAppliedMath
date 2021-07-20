%% Clear memory
%clc;
clear;
close all;

%Tipo de entrada (0) txt (1) CAD
flag_input = 1;
if (flag_input)
    %% Test model
    nFaces = 4;                    % number of faces from geometry
    nEl = 16;                     %numer of elements
    nO = 1;                        % nodes per element
    nNode = nO*nEl;              %number of nodes
    nodeInt = 2;                   %number of int nodes
    
    k = zeros(1,nEl);             % Conductivity
    len = zeros(1,nEl);           % array of elements lengths [m]
    q = ones(2,nFaces);             % array of prescibed contourn values [1] for unknown %mudar para variavel booleana
    nodesCor = zeros(2,nNode);      % array for nodes coordinates [x] and [y]
    elBC = zeros(3,nEl);        % face and boundary condition for each element
    nelCord = zeros(4,nEl);        % Coord of begin and start for each element
    intCord = zeros(2,nodeInt);     % coordinates for each int nod
    
    q(2,1) = 0;                 % Boundary condition per face ...
    q(1,2) = 200;               %[1] temperature [2] gradient
    q(2,3) = 0;
    q(1,4) = 100;
    
    for i = 1: nNode
        if i > 0 && i < 5
            nodesCor(1,i) = 0.125 + 0.25*(i-1);
            nodesCor(2,i) = 0;
        elseif i > 4 && i < 9
            nodesCor(1,i) = 1;
            nodesCor(2,i) = 0.125 + 0.25*(i-4-1);
        elseif i > 8 && i < 13
            nodesCor(1,i) = 1 - (0.125 + 0.25*(i-8-1));
            nodesCor(2,i) = 1;
        else
            nodesCor(1,i) = 0;
            nodesCor(2,i) = 1 -(0.125 + 0.25*(i-12-1));
        end
    end
    
    for i = 1:nEl
        k(i) = 1;
        len(i)= 0.25;
        if i > 0 && i < 5
            elBC(1,i) = 1;
            elBC(2,i) = q(1,1);
            elBC(3,i) = q(2,1);
            nelCord(1,i) = nodesCor(1,i) - len(i)/2;
            nelCord(2,i) = 0;
            nelCord(3,i) = nodesCor(1,i) + len(i)/2;
            nelCord(4,i) = 0;
        elseif i > 4 && i < 9
            elBC(1,i) = 2;
            elBC(2,i) = q(1,2);
            elBC(3,i) = q(2,2);
            nelCord(1,i) = 1;
            nelCord(2,i) = nodesCor(2,i) - len(i)/2;
            nelCord(3,i) = 1;
            nelCord(4,i) = nodesCor(2,i) + len(i)/2;
        elseif i > 8 && i < 13
            elBC(1,i) = 3;
            elBC(2,i) = q(1,3);
            elBC(3,i) = q(2,3);
            nelCord(1,i) = nodesCor(1,i) + len(i)/2;
            nelCord(2,i) = 1;
            nelCord(3,i) = nodesCor(1,i) - len(i)/2;
            nelCord(4,i) = 1;
        else
            elBC(1,i) = 4;
            elBC(2,i) = q(1,4);
            elBC(3,i) = q(2,4);
            nelCord(1,i) = 0;
            nelCord(2,i) = nodesCor(2,i) + len(i)/2;
            nelCord(3,i) = 0;
            nelCord(4,i) = nodesCor(2,i) - len(i)/2;
        end
        
        intCord(1,1) = 1;
        intCord(2,1) = 0.5;
        intCord(1,2) = 0;
        intCord(2,2) = 0.5;
        
    end
else
    points = xlsread('Input','Planilha1');
    nEl = points(1,2);
    pointsCor = points(2:end,1:2)';
    typeBC =points(2:end,4);
    BC = points(2:end,5);
    elBC = zeros(3,nEl);
    nelCord = zeros(4,nEl);  
    k = zeros(1,nEl); 
    len = zeros(1,nEl);   
    nFaces = nEl;
    q = ones(2,nFaces);  
    nO = 1;
    nodeInt = 1;
    nNode = nEl;
    intCord = [0.5;0.5];
    for i = 1: nEl
        nelCord(1:2,i) = pointsCor(:,i);
        if (i==nEl)
            nelCord(3:4,i) = pointsCor(:,1);
            r = (nelCord(3:4,i) - nelCord(1:2,i));
            len(i) = (dot(r,r))^0.5;
            continue;
        end
        nelCord(3:4,i) = pointsCor(:,i+1);
        r = (nelCord(3:4,i) - nelCord(1:2,i));
        len(i) = (dot(r,r))^0.5;
    end
    for i = 1:nFaces
        q(typeBC(i),i) = BC(i);
        elBC(:,i) = [i;q(:,i)];
    end
    nodesCor = ([nelCord(1:2,:)] + [nelCord(3:4,:)])/2;
    
end
%% Solve model iteratively and print each stage
potential = PotentialSolver(nFaces,nO,nEl,nodesCor,k,len,elBC,nelCord, ...
    nodeInt, intCord);
potential.DrawMesh();
potential.Solve();
potential.Solveint();


