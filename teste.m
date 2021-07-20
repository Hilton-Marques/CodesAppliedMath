clc;
clear all;
close all;
n = 50;
Xi = 10 * rand(2,n) - 5;
Xf = 10 * rand(2,n) - 5;
for i=1:n
    h=line(Xi(:,i),Xf(:,i));
    axis([-5 5 -5 5])
    pause(.1)
    delete(h)
end