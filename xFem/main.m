% Clear workspace
clear
close(findall(0,'Type','figure'));
clc;
n = 129; 
a = -1;
b = -1;
w = 2;
h = 2;
x = linspace(a,a+w,n);
y = linspace(b,b + h,n);
center = [a + w/2; b + h/2];
[X,Y] = meshgrid(x,y);
malha((n-1)^2) = Quadrado;
count = 1;
for i = 1:n-1
    for j = 1: n-1
        p1 = [X(i,j),Y(i,j)];
        p2 = [X(i,j+1),Y(i,j+1)];
        p3 = [X(i+1,j+1),Y(i+1,j+1)];
        p4 = [X(i+1,j),Y(i+1,j)];
        malha(count) = Quadrado([p1;p2;p3;p4],count);
        count = count + 1;
    end
end
nC = 8000;
t = linspace(0,2*pi,nC);
circle = center + 0.4*[cos(t);sin(t)];
%main loop
points = []; % pontos de interesse
points(end+1,:) = circle(:,1);
ids = []; % ids dos elementos
for i = 1:nC-1
    p1 = circle(:,i);
    p2 = circle(:,mod(i,nC-1)+1);
    id1 = [];
    id2 = [];
    for j = 1:length(malha)
        el = malha(j);
        if el.isInside(p1)
            id1 = el.id;
        end
        if el.isInside(p2)
            id2 = el.id;
        end
        if (~isempty(id1) && ~isempty(id2))
            break;
        end
    end
    if id1 == id2
        continue
    else
        heds = zeros(8,4);
        heds(1:4,:) = malha(id1).heds;
        heds(5:8,:) = malha(id2).heds;
        edgeCircle = [p1,p2];
        pt = intercession(heds,edgeCircle);
        points(end+1,1:2) = pt;
        ids(end+1) = id1;
        ids(end+1) = id2;
    end
end
%Plot
%points = round(points,6);
hold on
for i = 1:length(ids)
    el = malha(ids(i));
    el.plot;
end
plot(circle(1,:),circle(2,:),'Color','red');
plot(points(:,1),points(:,2),'o','Color','green');
fileID = fopen(strcat('points',num2str(n-1),'.txt'),'w');
fprintf(fileID,'%10s %10s\n','x','y');
for i = 1:length(points)
    p = points(i,:);
    fprintf(fileID,'%10.8f %10.8f\n',[p(1),p(2)]);
end
keyboard
function pt = intercession(heds,edgeCircle)
C = edgeCircle(:,1)';
D = edgeCircle(:,2)';
for i = 1:8
    hed = heds(i,:);
    A = hed(1:2) ;
    B = hed(3:4) ;
    if orient(C,D,A) > 0 && orient(C,D,B) > 0
        continue
    elseif orient(C,D,A) < 0 && orient(C,D,B) < 0
        continue
    elseif orient(A,B,C) > 0 && orient(A,B,D) > 0
        continue
    elseif orient(A,B,C) < 0 && orient(A,B,D) < 0
    end
    t = orient(C,D,A)/(orient(C,D,A) - orient(C,D,B));
    pt = (1 - t)*A+ t*B;
    return
end
end
function out = orient(A,B,C)
out = det([B-A;C-A]);
end
