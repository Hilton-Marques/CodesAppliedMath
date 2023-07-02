clear
close(findall(0,'Type','figure'));
clc;

N = 20;
u = linspace(0,2*pi,N+1);
v = linspace(-pi/2,pi/2,N+1);
[X,V] = meshgrid(u,v);
pts = zeros(N^2,3);
quad = Quad;
count = 1;
countQuad = 1;
hold on
view(30,30);
for j = 1:N
    for i = 1:N
        tetas = [u(i), u(i+1), u(i+1) , u(i)];
        fis = [v(j), v(j), v(j+1), v(j+1)];
        id = zeros(1,4);
        for k = 1: 4
            pts(count,:) = [cos(tetas(k))*cos(fis(k)), sin(tetas(k))*cos(fis(k)), ...
                sin(fis(k))];
            id(k) = count;
            count = count + 1;
        end
        quad.indices = id;
        
        quad.plot(pts);
    end
end
a = 1;

