clear
close(findall(0,'Type','figure'));
clc;

h = 1;
L = 3;

beam = Beam(h,L);
theta_f = pi/3;
len_final = (2*theta_f)*L/(2*sin(theta_f));
n = 100;
thetas = linspace(0,theta_f,n);
fig = figure ;
hold on
axis equal
axis off
set(gcf,'color','w');
axis(beam.getBB());
im = {};
for theta = thetas
    beam.updateCurvature(theta);
    beam.plot();
    pause(0.01);
    frame = getframe(fig);
    im{end+1} = frame2im(frame);
    beam.clean();
end
beam.plot();
frame = getframe(fig);
im{end+1} = frame2im(frame);
filename = 'beam.gif';
for idx = 1:length(im)
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1.0);
    elseif (idx > 1) && (idx < length(im))
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.04);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',3);
    end
end
