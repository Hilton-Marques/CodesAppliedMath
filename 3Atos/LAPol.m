clc;
clear all;
close all;
x = linspace(-pi,pi,100);
f1 = linspace(1,1,100);
f2 = x;
f3 = x.^2;
f4 = x.^3;
f5 = x.^4;
f6 = x.^5;
f7 = sin(x);
f8 = 0.987862*x - 0.155271*x.^3 + 0.00564312*x.^5;
plot(x,f1)
axis([-pi pi    -1 1 ])
pause(1)
hold on
plot(x,f2)
pause(1)
hold on
plot(x,f3)
pause(1)
hold on
plot(x,f4)
pause(1)
hold on
plot(x,f5)
pause(1)
hold on
plot(x,f6)
pause(5)
hold on 
plot(x,f7,'-o','MarkerIndices',1:5:length(f7))
hold on
legend('sin(x)','x','x^2','x^3','x^4','x^5','sin(x)')
pause(1)
n = 50;
for i=1:n
    polt = rand*f1 + rand*f2 + rand*f3 + rand*f4 + rand*f5 + rand*f6;
    h = plot(x,polt,'-c*','MarkerIndices',1:5:length(polt));
    pause(1)
    delete(h)
end
pause(1)
plot(x,f8,'-c*','MarkerIndices',1:5:length(f8),'color',[0 1 1])
legend('sin(x)','aprox','x^2','x^3','x^4','x^5','sin(x)')