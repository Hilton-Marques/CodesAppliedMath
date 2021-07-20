

x = 0:0.1:1;
y = 0:0.1:1;
z = 0:0.1:1;
dx = diff(x);
dy = diff(y);
dL = hypot(dx,dy);
[X,Y] = meshgrid([-1 1]);
cla
plot3(x,y,z,'-')
axis equal
hold on
for i = 1:length(x)-1 -4000         % too much rectangles
    X1 = X/2;                       % rectangle width is 1
    Y1 = Y*dL(i)/2;                 % rectangle height is dL
        % rotation matrix
    R = [dx(i) -dy(i); dy(i) dx(i)]/dL(i);
    V = R*[Y1(:) X1(:)]';           
    X1 = reshape(V(1,:),[2 2])+x(i)+dx(i)/2;
    Y1 = reshape(V(2,:),[2 2])+y(i)+dy(i)/2;
    Z = [z(i) z(i); z(i+1) z(i+1)];
    surf(X1,Y1,Z)
     pause(0.5)
end
hold off