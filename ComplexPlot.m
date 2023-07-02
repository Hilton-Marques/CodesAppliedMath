clc;
clear all;
close all;


%% Plot ln(1-z)vs Series Aproach map
hold on
n = 10;
z = complex(0,0.9);
eV = ln_1_zRender(z,n);
err = abs((eV(n) - log(1-z)))
hold off;
keyboard;

%% Maps
n = 250;   % number of steps
nVoltas = 1;  %number of voltas
dTheta = 0.1;    %step
%Initialization
i = 2;
zi = zeros(1,n);  %Domain
w = zeros(1,n);   %Image
iVoltasThroughOrigin = 0;
k = 3;          % flag for number of voltas
z0 = complex(1,1);  %initial value
z_aux0 = z0/8;     % auxiliar value to create path
%w(1) = ln_z(z0,0);  % map
w(1) = z_powTo(z0,1/3,0);  % map
zi(1) = z0;
z_aux = z_aux0;
%Plot initial vectors
figure('Name','Branch Points','NumberTitle','off','Position', [10 10 900 600]);
subplot(2,1,1);
title('Range z^(^1^/^3^)');
%Image
%Get BoundBox
w_max = z_powTo(z0,1/3,nVoltas);
%w_max = ln_z(z0,(n/(2*pi/dTheta)));
x_max = abs(real(w_max)) + 1.0;
x_min = -abs(real(w_max)) - 1.0;
y_max =  abs(imag(w_max)) + 1.0;
y_min =  -abs(imag(w_max)) - 1.0;
hold on
axis([x_min x_max y_min y_max])
vector2D([real(w(1)); imag(w(1))],[0 1 0]);
subplot(2,1,2);
title('Domain z^(^1^/^3^)');
%Domain
z_max = (n/(2*pi/dTheta))*z0;
x_max = abs(real(z_max));
x_min = -abs(real(z_max));
y_max =  abs(imag(z_max));
y_min =  -abs(imag(z_max));
hold on
axis([x_min x_max y_min y_max])
vector2D([real(zi(1)); imag(zi(1))],[0 1 0]);
pause(2);
while( i <= n)
    zi(i) = exp(complex(0,i*dTheta))*z_aux + (z0 - z_aux);
    z_relative = conj(zi(i))*z0;
    if (sign(angle(zi(i))) == 1 && (angle(zi(i))) <= dTheta )
        iVoltasThroughOrigin = iVoltasThroughOrigin + 1;
    end
    %w(i) = ln_z(zi(i),iVoltasThroughOrigin);
    w(i) = z_powTo(zi(i),1/3,iVoltasThroughOrigin);
    subplot(2,1,1);
    %Test if the zi pass to initial value
    if ( abs((real(z_relative) - abs(z0)^2)) <= 0.01 ...
            && sign(imag(z_relative)) == 1 )
        z_aux = k*z_aux0;
        vector2D([real(w(i)); imag(w((i)))],[0 1 0]);
        
        if (k < nVoltas)
            iVoltasThroughOrigin = 0;
            k = k + 1;

        end
    end
    vector2D([real(w(i));imag(w(i))],(i/n)*[0 0 1],[real(w(i-1));imag(w(i-1))])
    subplot(2,1,2);
    vector2D([real(zi(i));imag(zi(i))],(i/n)*[1 0 0],[real(zi(i-1));imag(zi(i-1))])
    i = i + 1;
    pause(0.05);
end
keyboard;

%% Plot (1+z)^n vs SeriesAproach map
hold on
n = 1/3;
m = 20;
z = complex(0,2);
eV = binomialTheoremRender(z,n,m);
w = (1 + z )^(1/3);
vector2D([real(w);imag(w)],(i/n)*[0 0 1],[real(w);imag(w)])
keyboard;

%% Plot exp(z)vs Aproach map
hold on
n = 5;
z = complex(0,-1);
eV = expRender(z,n);
z = complex(1,-1);
pause(1);
eV = expRender(z,n);
err = (eV(n) - exp(z));
hold off;
keyboard


x = -1:0.2:1;
[X,Y]=meshgrid(x);
Of = complex(X,Y);



%Spiral by a continuous complex time function
n = 10;    % number discrete intervals
tf = 4*pi; % final time
ti = 0;    % begin time
uo = 1; u = uo;  % initial positin
%ContinuousSpiral(ti,tf,n,uo);

%Spiral create by one complex number with translation conjugation
n = 20;              % number interations
source = 1.1*exp(complex(0,pi/10)); %transformation
uo = (3*exp(complex(0,pi/6)));      % initial pos
change = complex(1,-1);
%DiscreteSpiral(source,change,uo,n);

%Plane motion
%Final obj
ObjF = zeros(2,1);
I = zeros(1,1);
%Obj
c=0.2;
d = 1.5;
M = [ -1,(-1 + c),(-1 +2*c),(-1 +3*c),(-1 +4*c),(-1 +5*c),(-1 +6*c),(-1 +7*c) ...
    (-1 +8*c),(-1 +9*c),(-1 +10*c),(-1 +10*c),1.5,1,1,-1,-1 ; ...
    2,(2-c), 2,(2-c),2,(2-c),2,(2-c),2,(2-c),2,(2 - d),-1,-1,-2,-2,2];
O1 = complex(M(1,:),M(2,:));
hold on
O1 = Translation(false,O1,complex(1,3));
Conjugate(true,O1,pi/4,complex(1,-1));
keyboard;
%T2 = Rotation(0,(Translation(0,(Scale(0,Inversion(0,O1),4)),[-4;-4])),pi);
%Transformation(T2,O1);
t = 0:0.1:2*pi;
O2 = 1*exp(complex(0,t));

Id(1) = size(O2,2);
%Id(2) = size(O2,2);
O = [O2];
change = complex(2,0);
%DrawObject(O);
T = @(z) 2*z + change;
%T(O);
%Transformation(T,O,Id)



function triangle(xf,c,xi)
x = xf - xi;
L = vecnorm(x);
factorReduc = 1/5;
h = factorReduc*0.25;
b = factorReduc*0.25;
p_x = [-x(2,1); x(1,1)];
p_x = p_x/norm(p_x);
p1 = xf + b/2*p_x;
p2 = xf - b/2*p_x;
p3 = xf + h*(x/norm(x));

X = [p1(1,1), p2(1,1) , p3(1,1)];
Y = [p1(2,1), p2(2,1) , p3(2,1)];
Thandle = fill(X, Y, c);
end

function vector2D(xf,c,xi)
if (nargin < 2)
    c = [1 0 1];
end
if (nargin < 3)
    xi = zeros(2,1);
end
X = [xi(1,1), xf(1,1)];
Y = [xi(2,1), xf(2,1)];
line( X, Y, 'Color', c,'Linewidth',2);
triangle(xf,c,xi);
end

function DiscreteSpiral(source,change,uo,n)
un = uo + change;
power = 0:n+1;
T = source*ones(1,n);
P = diag(T.^(power'));
u = P'.*(uo * ones(1,n) );
u = P'.*((uo)*ones(1,n)) + change*ones(1,n);
plot(u);
hold on
vector2D([real(un);imag(un)], [0 1 1], ...
    [real(change);imag(change)]);
hold on
vector2D([real(change);imag(change)], [1 0 1]);
grid on
axis([-10 10 -10 10]);
drawnow;
end

function ContinuousSpiral(ti,tf,n,uo)
t = linspace(ti,tf,n);
dt = t(2) - t(1);
u = uo;
for tn = 1:ceil(tf/dt)
    Ri = ti:0.1:tn*dt;
    u = Ri.*exp(complex(0,Ri));
    clf;
    plot(u,'LineWidth',4);
    hold on
    %u = (-Ri).*exp(complex(0,Ri))*uo;
    %plot(u,'LineWidth',4);
    axis equal;
    axis([-10 10 -10 10]);
    grid on
    drawnow;
    %pause(1);
end
end
function Bbox = GetBoundBox(O)
Bbox = zeros(4,1);
border = 2;
X = real(O);
Y = imag(O);
xmax = max(max(X)) + border;
ymax = max(max(Y)) + border;
xmin = -xmax;
ymin = - ymax;
Bbox(1,1) = xmin;
Bbox(2,1) = xmax;
Bbox(3,1) = ymin;
Bbox(4,1) = ymax;
end
function DrawObject(O)
plot(O)
Bb = GetBoundBox(O);
axis(Bb);
grid on;
drawnow;
end
function Transformation(T,O,Id)
[m,n] = size(O);
if isa(T,'function_handle')
    F = T(O);
else
    F = T;
end
if nargin < 3
    Id = [n];
end
I = complex(ones(m,n),zeros(m,n));
nFrames = 10;
S = abs(F)./abs(O);
R = exp(complex(0,angle(F)) - complex(0,angle(O)));
Bb = GetBoundBox(F);
%Bb = GetBoundBox(O);
Bb = [-5 5 -5 5];
plot(O ,'Color','red');
axis(Bb);
grid on;
pause(1);
for i = 1:nFrames
    Si = I + (i/nFrames)*(S - I);
    Ri =  R.^(i/nFrames) ;
    Fi = (Ri.*Si).*O;
    %Fi = I.*O + (i/nFrames)*(F - I.*O);
    k = 1;
    l = 1;
    for j = 1:size(Id,2)
        l = l + Id(j)-1;
        h = plot(Fi(1:end,k:l),'Color','blue');
        k = k + Id(j);
    end
    drawnow;
    pause(0.3);
    delete(h);
end
plot(F,'Color','blue');
end
function Circle(color)
if nargin < 1
    color = [1 0 0];
end
t = linspace(0,2*pi,50);
c = exp(complex(0,t));
plot(c,'Color',color);
end

function T = Translation(flag,O,offset)
T = O + offset;
vector2D([real(offset);imag(offset)],[0 1 1]);
if flag == false
    return
end
f = @(z) z + offset;
Transformation(f,O);
end

function T = Rotation(flag,O,angle,center)
if nargin < 4
    center = complex(0,0);
end
R = exp(complex(0,angle));
T = (R*O) + center;
vector2D([real(center);imag(center)],[0 1 1]);
if flag == false
    return
end
f = @(z) R*z + center;
Transformation(f,O);

end
function T = Scale(flag,O,factor)
T = factor*O;
if flag == false
    return
end
f = @(z) factor*z;
Transformation(f,O);
end
function T = Inversion(flag,O)
f = @(z) 1./z;
T = f(O);
if flag == false
    return
end
Transformation(f,O);
end
function T = Conjugate(flag,O,angle,b)
a = 1;
if nargin > 2
    a = exp(complex(0,angle));
end
if nargin < 4
    b = 0;
    off = 0;
else
    off = b - (a^2)*(abs(b))^2/b;
end
%Plot line reflection
scaleFactor = 10;
xf = b + scaleFactor*a;
xi = b - scaleFactor*a;
grid on
line([real(xi) real(xf)],[imag(xi) imag(xf)],'Color','blue','LineStyle','--');
vector2D([real(b);imag(b)],[0.5 0.5 0.5]);

f = @(z) (a^2)*(abs(z)).^2 ./z + off;
T = f(O);
if flag == false
    return
end
Transformation(f,O);
end

function expRender = expRender(z,n)
axis([-3 3 -3 3]);
expRender = zeros(n,1);
expRender(1) = 1;
vector2D([real(expRender(1));imag(expRender(1))],[0 0 1]);
pause(0.1);
i = 2;
while( i <= n)
    expRender(i) = expRender(i-1) + (z^(i-1)) / factorial((i-1));
    vector2D([real(expRender(i));imag(expRender(i))],[0 0 1],[real(expRender(i-1));imag(expRender(i-1))])
    i = i + 1;
    pause(0.5);
end
vector2D([real(exp(z));imag(exp(z))],[1 0 1]);
pause(0.1);
end
function bT = binomialTheoremRender(z,n,m)
axis([-0.5 1.5 -0.5 1]);
bT = zeros(m,1);
bT(1) = 1;
vector2D([real(bT(1));imag(bT(1))],[0 0 1]);
pause(0.1);
i = 2;
while( i <= m)
    bT(i) = bT(i-1) + (fact(n,i-1)*z^(i-1))/factorial(i-1) ;
    vector2D([real(bT(i));imag(bT(i))],[0 0 1],[real(bT(i-1));imag(bT(i-1))])
    i = i + 1;
    pause(0.5);
end
vector2D([real((1+z)^n);imag((1+z)^n)],[1 0 1]);
pause(0.1);
end
function out = fact(n,m)
fac(1) = n;
i = 2;
if (m == 1)
    out = n;
    return;
end
while (i <= m)
    fac(i) = fac(i-1)*(n-1);
    i = i+1;
end
out = fac(m);
end
function lnSeries = ln_1_zRender(z,n)
Circle([0 0 0]);
axis([-1.3 1.3 -1.3 1.3]);
lnSeries(1) = -(z^(1)/1);
vector2D([real(lnSeries(1));imag(lnSeries(1))],[0 0 1])
pause(0.5)
i = 2;
while( i <= n)
    lnSeries(i) = lnSeries(i-1) - (z^(i)/i);
    vector2D([real(lnSeries(i));imag(lnSeries(i))],[0 0 1],[real(lnSeries(i-1));imag(lnSeries(i-1))])
    i = i + 1;
    pause(0.5);
end
vector2D([real(log(1-z));imag(log(1-z))],[1 0 1]);
pause(0.1);
end
function w = z_powTo(z,pow,nVoltas)
R = abs(z);
w_R = (R^pow);
if (angle(z) > 0)
    theta = angle(z) + nVoltas*2*pi;
    w_theta = theta*pow;
    w = w_R*exp(complex(0,w_theta));
else
    theta = angle(z) + 2*pi + nVoltas*2*pi;
    w_theta = theta*pow;
    w = w_R*exp(complex(0,w_theta));
end
end



function w = ln_z(z,nVoltas)
R = abs(z);
a = log(R);
if (angle(z) > 0)
    b = angle(z) + nVoltas*2*pi;
    w = complex(a,b);
else
    b = angle(z) + 2*pi + nVoltas*2*pi;
    w = complex(a,b);
end
end

