clear
close(findall(0,'Type','figure'));
clc;

t = linspace(0,2*pi,100);
curve = 2*abs(cos(2*t)).*[cos(t); sin(t)];
hold on
plot(curve(1,:), curve(2,:),'color','red');
plot(0,1,'x','MarkerSize',10);
plot(0,-1,'x','MarkerSize',10);

keyboard;
curve = (0.25)*[16*sin(t).^3 ; 13*cos(t) - 5*cos(2*t) - 2*cos(3*t) - cos(4*t)];
circle = [cos(t); sin(t)];
hold on
curve = curve ./ (1 + 1*(vecnorm(curve,1) - 1));
%quiver(0*curve(1,:),0*curve(1,:),curve(1,:),curve(2,:));
plot(curve(1,:), curve(2,:),'color','red');
s = linspace(0,1,5);
for i = 1:5
    temp = curve ./ (1 + s(i)*(vecnorm(curve,1) - 1));
    temp = (curve + s(i)*(circle - curve));
    %quiver(0*temp(1,:),0*temp(1,:),temp(1,:),temp(2,:));
    plot(temp(1,:), temp(2,:),'--');
end
plot(temp(1,:), temp(2,:),'color','blue')
