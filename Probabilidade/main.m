clc;
clear all;
close all;

figure
hold on
grid off
axis off
axis equal
xlim([-0.2,1.2]);
ylim([-0.2,1.2]);

q4();
keyboard
b1 = Box([0,0],[1,1],'S','red');
b2 = Box([0.1,0.1],[0.6,0.7],'A','blue');
b2.Divide();
e = Elipse(b2.getCenter(),b2.getD()/5,b2.getD()/5.5);
e.show();

function q1()
b1 = Box([0,0],[1,1],'','red');
b2 = Box([0,0],[0.6,1],'$x_{t}{$','blue');
b3 = Box([0,0],[1,0.68],'$z_{1:t-1},u_{1:t}$','green');
%b4 = Box([0,1],[1,0.61],'$z_{t}$','magenta');
exportgraphics(gca,'prob_2.png','Resolution',1000);
end
function q2()
b1 = Box([0,0],[1,1],'','red');
b2 = Box([0,0],[0.6,1],'$x_{t}$','blue');
b3 = Box([0,0],[1,0.66],'$z_{1:t-1},u_{1:t}$','green');
b1.Divide('$x_{t-1}^{');

exportgraphics(gca,'prob_3.png','Resolution',1000);

end
function q3()
fac = 0.1;
color =  [0    102/255   51/255];
%b1 = Box([0,0],[0.6,0.66],'',color);
%exportgraphics(gca,'prob_6.png','Resolution',1000);

b1 = Box([0,0],[0.5,0.5],'',color);
b2 = Box([0,0.5]+[0,fac],[0.5,0.66]+[0,fac],'',color);
b3 = Box([0.5,0.5]+[fac,fac],[0.6,0.66]+[fac,fac],'',color);
b4 = Box([0.5,0.0]+[fac,0.0],[0.6,0.5]+[fac,0.0],'',color);

exportgraphics(gca,'prob_7.png','Resolution',1000);

end
function q4()
fac = 0.1;
color =  [0    122/255   71/255];
color2 = [122/255, 224/255,0];
%b1 = Box([0,0],[0.6,0.66],'',color);
%b1 = Box([0.6,0],[1.0,0.66],'',color2);
%exportgraphics(gca,'prob_5.png','Resolution',1000);


b1 = Box([0,0],[0.5,0.5],'',color);
b2 = Box([0,0.5]+[0,fac],[0.5,0.66]+[0,fac],'',color);
b3 = Box([0.5,0.5]+[fac,fac],[0.6,0.66]+[fac,fac],'',color);
b3 = Box([0.6,0.5]+[fac,fac],[1.0,0.66]+[fac,fac],'',color2);

b4 = Box([0.5,0.0]+[fac,0.0],[0.6,0.5]+[fac,0.0],'',color);
b4 = Box([0.6,0.0]+[fac,0.0],[1.0,0.5]+[fac,0.0],'',color2);
exportgraphics(gca,'prob_8.png','Resolution',1000);

end