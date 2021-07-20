clc;
clear all;
close all;
L = 10;

for i = 1:5
x = linspace(-L,L,101);
F = f(x,L,i);
clf;
plot(x,F);
hold;
axis([-L L 0 1]);
drawnow;
pause(0.1);

end
function F = f(w,L,i)
F = 0.25*(1 + cos((pi/L)*w)).*(1+cos(w + i));
end