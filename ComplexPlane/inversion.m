clear
close(findall(0,'Type','figure'));
clc;

Mobius()
keyboard
euler = quat2eul([0.00666,0.0003784,0.05796,0.998296])
matrix = eul2rotm(euler)
x = matrix(:,1)
y = matrix(:,2)
goal = [4,4,0];
cur = [0.01423,0.56834,0]
x_coord = dot(x,goal - cur)
keyboard
global i 
i = complex(0,1);

%Invert semi plane
x = linspace(-10,10,10);
y = linspace(0,10,20);
[X,Y] = meshgrid(x,y);
z = X + complex(0,1)*Y;
z = linspace(0,100,100);
mobius = (z - i ) ./ (z + i);
figure
hold on 
plot_z(z);
hold off
figure
hold on 
plot_z(mobius);
teta = linspace(0,2*pi,61);
circle = exp(teta*i);
plot(real(circle), imag(circle),'color','black')
%plot_z(circle)
hold off
keyboard
teta = linspace(0,2*pi,31);
circle = exp(teta*i);
fig = figure;
hold on
axis ([-2,2,-2,2])
plot(0,0,'x','MarkerSize',10)
%plot_nice_circle()
plot_ortho_circles()
keyboard
plot(real(circle), imag(circle))
a = complex(0.5,0.5);
a_til = 1/conj(a);
plot_z(a)
plot_z(a_til)
plot_z(a*a_til)
a_reta = complex(-0.72,0.68);
plot_line(a_reta);
z = line_z(a_reta);
z_til = 1./conj(z);
plot(real(z_til),imag(z_til));
z = complex(-0.3,0.6);
z_ref = reflection(complex(0,0),z);
plot_z(z);
plot_z(z_ref);
keyboard;

function plot_z(z)
    plot(real(z),imag(z),'o');
end
function plot_line(a)
t = linspace(-1,1,20);
z = a + complex(0,1) * (1 / conj(a)) * t;
plot(real(z),imag(z));
end
function z =line_z(a)
t = linspace(-10,10);
z = a + complex(0,1) * (1 / conj(a)) * t;
end
function z_out = reflection(a,z)
    z_out = z + a;
    return
    teta = angle(i*a);
    z_out = exp(complex(0,teta))*conj(exp(complex(0,-teta))*z) + a;
end
function plot_nice_circle()
    teta = linspace(0,2*pi,31);
    circle = exp(complex(0,1)*teta);
    plot(real(circle),imag(circle));
    t = linspace(0,1,30);
    for ti = t
        h = [];
        w = [];
        center = -ti*complex(0,1);
        %w = quiver(0,0,real(center),imag(center),'color','cyan');
    for z = circle
        sinal = 1;
        if real(z) < 0 
            sinal = -1;
        end        
        z_til = z - center;
        mod_z = abs(z_til);
        fac = sqrt(4 - mod_z^2)/mod_z
        y = z + fac*(sinal)*complex(0,1)*z_til;    
        l = [z,y];        
        h(end+1) = plot(real(l),imag(l));        
    end
        pause(0.2);
        delete(w);
        delete(h);
        
    end
end
function plot_nice_inverse()
    teta = linspace(0,2*pi,31);
    circle = exp(complex(0,1)*teta);
    plot(real(circle),imag(circle));
    fixed_z = complex(0,1);
    t = linspace(0,1,30);
    for ti = t
        center = complex(0,ti);
        R = abs((fixed_z - center))*0.5;
        center_circle = center + (fixed_z - center)*0.5;
        new_circle = R*circle + center_circle;
        z_til = 1 ./ conj(new_circle);
        h = plot(real(new_circle), imag(new_circle));
        w = plot(real(z_til),imag(z_til));
        pause(1);     
        delete(h);
        delete(w);
    end
end
function plot_ortho_circles()
    teta = linspace(0,2*pi,31);
    circle = exp(complex(0,1)*teta);
    plot(real(circle),imag(circle));
    for z_fixed = circle
        dir = z_fixed + 3*complex(0,1)*z_fixed;
        dir_opp = z_fixed - 3*complex(0,1)*z_fixed;
        line = [dir_opp,dir];
        plot(real(line),imag(line));
        z_teste = complex(0.5,0.5);
        z_teste_til = 1 / conj(z_teste);
        plot_z(z_teste);
        plot_z(z_teste_til);
        corda = z_teste_til - z_teste;
        corda_len = abs(corda);
        mid = z_teste + corda * 0.5;
        other_dir = mid + 2*complex(0,1)*mid;
        other_dir_opp = mid - 2*complex(0,1)*mid;
        line_2 = [other_dir_opp,other_dir];
        %plot(real(line_2),imag(line_2));
        center_d = inter(dir_opp,dir,other_dir_opp,other_dir);
        radius = abs(center_d - z_teste);
        new_circle = radius*circle + center_d;
        h = plot(new_circle);
        pause(1);
        %delete(h);
    end
end
function center = inter(z1,z2,z3,z4)
    A = [real(z1),imag(z1)];
    B = [real(z2),imag(z2)];
    C = [real(z3),imag(z3)];
    D = [real(z4),imag(z4)];
    t = orient(C,D,A)/(orient(C,D,A) - orient(C,D,B));
    pt = (1 - t)*A+ t*B;
    center = complex(pt(1),pt(2));
end
function out = orient(A,B,C)
out = det([B-A;C-A]);
end
function out = unit_circle()
teta = linspace(0,2*pi,20);
rot = exp(complex(0,teta));
radius = linspace(0.1,1,5);
out = kron(rot,radius);
end
function Mobius()
%mobius
q = exp(complex(1,0));
r = exp(complex(0,(2*pi/3)));
s = exp(complex(0,(4*pi/3)));
z = unit_circle();
x = linspace(-100,100,100);
y = linspace(-1,-100,100);
[X,Y] = meshgrid(x,y);
z = complex(X,Y);
z = z(:);
z = unit_circle();
z = 0.99*z
mobius = complex(0,1) * (1 - z) ./ (z + 1) - complex(0,1)
mobius = z.^2
%mobius = ((z - q)* (r - s)) ./ ((z-s) .* (r-q));
mobius_1 = (complex(0,1) - z) ./  (z + complex(0,1));
mobius_2 = (z - complex(0,1)) ./  (z + complex(0,1));
%mobius = [mobius_1;mobius_2];
%mobius = (1-i).*(z-i)./(z-1)
figure
%axis([-10,10,-1,1])
hold on 
plot_z(z);
hold off
figure
hold on 
axis([-1,1,-1,1])
plot_z(mobius);
hold off
end
function out = mobius_inv(z)
out = ((-1)+(-1).*i+((-1)+i).*z).^(-1).*(1+i+((-1)+i).*z);
end