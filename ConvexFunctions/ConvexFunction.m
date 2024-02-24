clc;
clear all;
close all;


%Default configurations
fig = figure;
hold on
angleInit = 90;
view(angleInit,0);
axis ([-10,10,-10,10]);
set(gca,'visible','off');
set(gcf,'color','white');
count = 1;
frame = getframe(fig);
im{count} = frame2im(frame);

x = linspace(-10,10,100);
y = linspace(-10,10,100);
[X,Y] = meshgrid(x,y);
Z = X.^2 + Y.^2;
h = surf(X,Y,Z);

[fRot,count,im] = rot(angleInit,count,im,fig);

delete(h);
axis ([-10,10,-10,7]);
contour(X,Y,Z);
count = count + 1;
frame = getframe(fig);
im{count} = frame2im(frame);
%Plot x1
x1 = [-3.93939,-1.91919];
plot(x1(1),x1(2),'o','MarkerFaceColor','cyan');
count = count + 1;
frame = getframe(fig);
im{count} = frame2im(frame);
%Plot x2
x2 = [0.505051,-6.36364];
plot(x2(1),x2(2),'o','MarkerFaceColor','red');
count = count + 1;
frame = getframe(fig);
im{count} = frame2im(frame);
% u
u = x2 - x1;
quiver(x1(1),x1(2),u(1),u(2),'color','green');
count = count + 1;
frame = getframe(fig);
im{count} = frame2im(frame);
%Plot gradient 1
g1 = x1;
quiver(x1(1),x1(2),g1(1),g1(2),'color','cyan');
count = count + 1;
frame = getframe(fig);
im{count} = frame2im(frame);
%Plot gradient 2
g2 = 0.7*x2;
quiver(x2(1),x2(2),g2(1),g2(2),'color','red');
count = count + 1;
frame = getframe(fig);
im{count} = frame2im(frame);
%Translate
[count,im] = translate(g1,x1, [5,5],'cyan',im,count,fig);
[count,im] = translate(g2,x2, [5,5],'red',im,count,fig);
%Diff
diff = g2 - g1;
quiver(5,5,diff(1),diff(2),'color','blue');
count = count + 1;
frame = getframe(fig);
im{count} = frame2im(frame);

[count,im] = translate(diff,[5,5],x1,'blue',im,count,fig);

filename = 'convex_function.gif';
createGif(im,filename);

keyboard

function [fRot,count,im] = rot(initangle,count,im,fig)
count = count + 1;
frame = getframe(fig);
im{count} = frame2im(frame);
for i = 1:1:90
    fRot = initangle + i;
    view(fRot,i);
    pause(0.05);
    count = count + 1;
    frame = getframe(fig);
    im{count} = frame2im(frame);
end
end
function [count,im] = translate(v,come, fim,color,im,count,fig)
t = linspace(0,1,10);
for i = 1:10
   x_i = come + (fim - come) * t(i);
   h = quiver(x_i(1),x_i(2),v(1),v(2),'color',color);
   count = count + 1;
   frame = getframe(fig);
   im{count} = frame2im(frame);
   delete(h);
   pause(1);
end
quiver(x_i(1),x_i(2),v(1),v(2),'color',color);
count = count + 1;
frame = getframe(fig);
im{count} = frame2im(frame);
end
function createGif(im,filename)
for idx = 1:length(im)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.5);
    elseif (idx > 1) && (idx < 92)
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
    elseif (idx == 92)
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
    elseif (idx == 93)
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',2.5);
    elseif (idx > 93) && (idx < length(im))
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.5);
    elseif (idx == length(im))
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',5);
    end
end
end