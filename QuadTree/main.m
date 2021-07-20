% Clear workspace
clear
close(findall(0,'Type','figure'));
clc
%%Get Boundary
%points = xlsread('Points','Planilha1');
n = 5000;
points = zeros(n,2);
for i = 1:n
    points(i,:) = [rand, rand];
end
draw_obj(points);


x_max = max(points(:,1));
x_min = min(points(:,1));
y_min = min(points(:,2));
y_max = max(points(:,2));
x_c = (x_max + x_min)/2;
y_c = (y_max + y_min)/2;
w = ceil(x_max - x_min)/2;
h = ceil(y_max - y_min)/2; 
l = max(w,h);
w = l;
h = l;
%% Initialize QtTree
capacity = 1;
boundary = Rectangle(x_c,y_c,w,h);
qt = QuadTree(boundary,capacity);
tic
for i = 1:size(points,1)
    pt = points(i,:);
    qt.insert(pt,i);
end
toc
qt.show();


function draw_obj(points)
plot(points(:,1),points(:,2),'o');
for i = 1:size(points,1)
 text(points(i,1),points(i,2),num2str(i),'HorizontalAlignment','left');
end
end
