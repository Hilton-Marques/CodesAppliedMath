clear
close(findall(0,'Type','figure'));
clc;

%rotationScene('rotation')
scale('scale');

function out = f(x,y)
%out = 1+10.*y+4.*x.*y+40.*x.^5.*y+20.*x.*y.^3;
%out = ((-11)+x.^2+y).^2+((-7)+x+y.^2).^2;
out = -((-0.15E1).*x+(x+(-1).*y).^2+0.25E1.*y+sin(x+y));
end
function M = rotation(teta)
M = [[cos(teta),-sin(teta)];[sin(teta),cos(teta)]];
end
function out = laplace(x,y)
out = 4 - 2*sin(x+y);
end
function rotationScene(type)
m = 30;
x = linspace(-10,10,m);
[X_1,Y_1] = meshgrid(x,x);
Z = f(X_1,Y_1);
rot_matrix = rotation(pi/4);
p = [X_1(:)';Y_1(:)'];
p = rot_matrix * p;
X_2 = reshape(p(1,:),m,m);
Y_2 = reshape(p(2,:),m,m);

%getBB
figure
hold on
view(30,30);
set(gca,'XColor', 'none','YColor','none','ZColor','none');
surf(X_1,Y_1,Z);
surf(X_1,Y_1,-450*ones(size(X_1)),'FaceColor',[0.5,0.5,0.5],'EdgeAlpha',0);
surf(X_2,Y_2,Z);
box = [get(gca,'xlim'),get(gca,'ylim')];
hold off

figure
hold on
view(30,30);
set(gca,'XColor', 'none','YColor','none','ZColor','none');
surf(X_1,Y_1,Z);
surf(X_1,Y_1,-450*ones(size(X_1)),'FaceColor',[0.5,0.5,0.5],'EdgeAlpha',0);
axis(box);
exportgraphics(gca,strcat('inital_pose',type,'.jpeg'),'Resolution',333)

hold off

figure
hold on
view(30,30);
set(gca,'XColor', 'none','YColor','none','ZColor','none');
%Z = f(X,Y_1);
surf(X_2,Y_2,Z);
surf(X_2,Y_2,-450*ones(size(X_2)),'FaceColor',[0.5,0.5,0.5],'EdgeAlpha',0);
axis(box);
exportgraphics(gca,strcat('rotated',type,'.jpeg'),'Resolution',333)

hold off

Z_2 = laplace(X_2,Y_2);

figure
hold on
h = surf(X_2,Y_2,Z_2);


figure
hold on
view(30,30);
set(gca,'XColor', 'none','YColor','none','ZColor','none');
Z_2 = laplace(X_2,Y_2);
h1 = surf(X_2,Y_2,Z);
h1.CData = (h1.CData - mean(h1.CData,'all'))/std2(h1.CData);
cdata = (h.CData - mean(h.CData,'all'))/std2(h.CData);
surf(X_2,Y_2,-450*ones(size(X_2)),cdata);
axis(box);
exportgraphics(gca,strcat('final',type,'.jpeg'),'Resolution',333)

hold off


figure
hold on
view(30,30);
set(gca,'XColor', 'none','YColor','none','ZColor','none');
h1 = surf(X_1,Y_1,Z);
h1.CData = (h1.CData - mean(h1.CData,'all'))/std2(h1.CData);
surf(X_1,Y_1,-450*ones(size(X_2)),cdata);
axis(box);
exportgraphics(gca,strcat('laplacian',type,'.jpeg'),'Resolution',333)
hold off

end

function scale(type)
m = 30;
x = linspace(-10,10,m);
[X_1,Y_1] = meshgrid(x,x);
Z = f(X_1,Y_1);
lam = 2;
X_2 = lam*X_1;
Y_2 = lam*Y_1;

%getBB
figure
hold on
view(30,30);
set(gca,'XColor', 'none','YColor','none','ZColor','none');
surf(X_1,Y_1,Z);
surf(X_1,Y_1,-450*ones(size(X_1)),'FaceColor',[0.5,0.5,0.5],'EdgeAlpha',0);
surf(X_2,Y_2,Z);
box = [get(gca,'xlim'),get(gca,'ylim')];
hold off

figure
hold on
view(30,30);
set(gca,'XColor', 'none','YColor','none','ZColor','none');
surf(X_1,Y_1,Z);
surf(X_1,Y_1,-450*ones(size(X_1)),'FaceColor',[0.5,0.5,0.5],'EdgeAlpha',0);
axis(box);
exportgraphics(gca,strcat('inital_pose',type,'.jpeg'),'Resolution',333)

hold off

figure
hold on
view(30,30);
set(gca,'XColor', 'none','YColor','none','ZColor','none');
%Z = f(X,Y_1);
surf(X_2,Y_2,Z);
surf(X_2,Y_2,-450*ones(size(X_2)),'FaceColor',[0.5,0.5,0.5],'EdgeAlpha',0);
axis(box);
exportgraphics(gca,strcat('rotated',type,'.jpeg'),'Resolution',333)

hold off

Z_2 = laplace(X_2,Y_2);

figure
hold on
h = surf(X_2,Y_2,Z_2);


figure
hold on
view(30,30);
set(gca,'XColor', 'none','YColor','none','ZColor','none');
Z_2 = laplace(X_2,Y_2);
h1 = surf(X_2,Y_2,Z);
h1.CData = (h1.CData - mean(h1.CData,'all'))/std2(h1.CData);
cdata = (h.CData - mean(h.CData,'all'))/std2(h.CData);
surf(X_2,Y_2,-450*ones(size(X_2)),cdata);
axis(box);
exportgraphics(gca,strcat('final',type,'.jpeg'),'Resolution',333)

hold off


figure
hold on
view(30,30);
set(gca,'XColor', 'none','YColor','none','ZColor','none');
h1 = surf(X_1,Y_1,Z);
h1.CData = (h1.CData - mean(h1.CData,'all'))/std2(h1.CData);
surf(X_1,Y_1,-450*ones(size(X_2)),cdata);
axis(box);
exportgraphics(gca,strcat('laplacian',type,'.jpeg'),'Resolution',333)
hold off
end