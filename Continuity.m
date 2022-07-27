clear
close(findall(0,'Type','figure'));
clc;


figure
hold on
axis([-2,2,-1,2]);
z = complex(2,2);
out = conj(z) / z;
line([-2,2],[0,0],'LineStyle','--');
quiver(0,0,real(z),imag(z),'color','blue');
quiver(0,0,real(out),imag(out),'color','red');
DrawAngle(angle(z))
DrawAngle(angle(out),z,'-3')
DrawConformal()
keyboard
figure
hold on
teta = linspace(0,2*pi);
radius = linspace(2,0.1,10);
for r = radius
    circle = r*exp(complex(0,teta));
    for z = circle
        v = 0.1*- z ;
        quiver(real(z),imag(z),real(v),imag(v),'color','black');        
    end
    out = conj(circle) ./ circle;
    plot(real(circle),imag(circle),'blue');
    plot(real(out),imag(out),'red','linewidth',5);
    pause(1);
end
function DrawAngle(phase_angle,init,string)
if nargin == 1
    init = complex(1,0);
    string = '\theta';
end
if nargin == 3
    string = strcat(string,'\theta');
end
radius = 0.5 * 1/abs(init);
init_angle = angle(init);
teta = linspace(init_angle, phase_angle);
final_pos = radius*exp(complex(0,phase_angle));
arc = radius*exp(complex(0,teta));
plot(real(arc),imag(arc),'color','black');
pos_text = arc(ceil(end/2));
text(real(pos_text),imag(pos_text),string,'HorizontalAlignment','left');
% draw arrow
sinal = sign(phase_angle-init_angle);
dir = complex(0,1)*final_pos;
normalized_dir = sinal*0.01*dir/abs(dir);
quiver(real(final_pos),imag(final_pos),real(normalized_dir),imag(normalized_dir),'color','black','LineWidth',3)
end
function DrawConformal()
figure
hold on
axis ([-2,2,-2,2])
z = complex(1,1);
plot(real(z),imag(z),'x')
trans_1 = complex(1,1);
trans_2 = complex(-2,1);
center_1 = z + trans_1;
center_2 = z + trans_2;
angle_1 = angle(-trans_1);
angle_2 = angle(-trans_2);
teta_1 = linspace(angle_1 - pi/3,angle_1 + pi/3);
teta_2 = linspace(angle_2 - pi/3,angle_2 + pi/3);
circle_1 = center_1 + abs(trans_1) * exp(complex(0,teta_1));
circle_2 = center_2 + abs(trans_2) * exp(complex(0,teta_2));
plot(real(circle_1),imag(circle_1),'color','blue');
plot(real(circle_2),imag(circle_2),'color','red');
% plot output
out_1 = conj(circle_1) ./ circle_1;
out_2 = conj(circle_2) ./ circle_2;
plot(real(out_1),imag(out_1),'color','blue');
plot(real(out_2),imag(out_2),'color','red');
end