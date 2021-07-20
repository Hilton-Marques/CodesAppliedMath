clear
close all
clc

ng = 5;
n = 8;
a1 = 1;
an = 2*(n)+1;
size = (a1 + an)*(n+1)/2;
M = zeros(size,1);
L = zeros(size,1);
% Momentos
yc = [0.625000, 0.375000, 0.1250000];
xs = [1.625, 2.375, 3.125];
xl = [1.825, 2.575, 3.625];
p1 = [1, 1, 0];
p2 = [1, 0, 0];
p3 = [0, 0, 0];
area = findArea(p1,p2,p3);
[X,Y,Wx,Wy] = triQuad(ng);
count = 1;
tic 
for n_ = 1:n + 1
    mTot = 2*(n_- 1) + 1; % total m for each n
    for m_ = 1:mTot
        for l = 1:ng
            for k = 1:ng
                xksi = (1 - X(k,l) - Y(k,l))*p1 + X(k,l)*p2 + Y(k,l)*p3;
                r = xksi - yc;
                R(n_,r,[]);
                M(count) = M(count) + R(n_,r,m_)*Wx(k)*Wy(l)*(area/0.5);
            end
        end
        count = count + 1;
    end
end
M;
S = Sb(n+1,xs-yc,[]);
int = sum(S.*M)
%M2L
r = xl - yc;
Sb(4,r,[]);
count = 1;
for n_ = 0:n 
    for m_ = -n_:n_
        for nc = 0: n
            for mc = -nc:nc
                nreal = n_ + nc;
                mcTot = 2*(nreal) + 1;
                nvirtual = nreal + 1; % conversion due to Matlab
                mreal = m_+mc; 
                mvirtual = mreal + (mcTot - 1)/2  + 1; % conversion due to Matlab
                L(count) = L(count) + ...
                    ((-1)^(n_))*Sb(nvirtual,r,mvirtual)*M(indexNM(nc,mc));
            end
        end
        count = count + 1;
    end
end
%Integral
r = (xs - xl);
RR = R(n+1,r,[]);
int = RR*L
toc
function plotTriangle(p1,p2,p3)
line([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],'Color','black');
line([p2(1),p3(1)],[p2(2),p3(2)],[p2(3),p3(3)],'Color','black');
line([p3(1),p1(1)],[p3(2),p1(2)],[p3(3),p1(3)],'Color','black');
xc = (p1 + p2 + p3)/3;
%text(xc(1),xc(2),num2str(i));
end
function out = R(n,y,m)
if n == 1
    out = 1;
    return;
end
r2 = dot(y,y);
%% Initialize Values
RValues(n) = SphericalHarmonics();
for i = 1:n
    RValues(i) = SphericalHarmonics(i);
end
%% Obtain N,N Values
RValues(1).addNNValue(1);
RValues(1).merge();
for i = 2:n
    Rnn = RValues(i-1).MValuesPos(end);
    value = -Rnn*(y(1) + complex(0,1)*y(2))/(2*(i-1));
    RValues(i).addNNValue(value);
end
%% Obtain N,M Values with M => 0
R10 = y(3)*RValues(1).MValuesPos(1);
RValues(2).addMValues(R10,1);
RValues(2).merge();
for i = 3:n
    obj = RValues(i);
    objAn2 = RValues(i-2);
    objAn1 = RValues(i-1);
    mTot = obj.m;   %Calcule o total de indices m
    mEnd = 1+(mTot-1)/2;  % Calcule quanto tem do centro até a borda
    for j = 1:mEnd-2
        MValue = (1/(((i-2) + (j-1) + 1)*(i-2 + 1 - (j-1))))*...
            ( (2*(i-2)+1)*y(3)*objAn1.MValuesPos(j) - ...
            r2*objAn2.MValuesPos(j) );
        obj.addMValues(MValue,j);
    end
    j = j + 1;
    MValue = (1/(((i-2) + (j-1) + 1)*(i-2 + 1 - (j-1))))*...
        (2*(i-2)+1)*y(3)*objAn1.MValuesPos(j);
    obj.addMValues(MValue,j);
    obj.merge();  % Find the negative values of m
end
%Quando olhei no mathematica os sinais estavam inversos
out = RValues(end).values;
if ~isempty(m)
    out = out(m);
else
    out = [RValues(:).values];
end
end
function [X,Y,Wx,Wy] = triQuad(N)
v = [0 0; 0 1; 1 0];
n=1:N;  nnk=2*n+1; A=[1/3 repmat(1,1,N)./(nnk.*(nnk+2))];
n=2:N; nnk=nnk(n); B1=2/9; nk=n+1; nnk2=nnk.*nnk;
B=4*(n.*nk).^2./(nnk2.*nnk2-nnk2); ab=[A' [2; B1; B']]; s=sqrt(ab(2:N,2));
[V,X]=eig(diag(ab(1:N,1),0)+diag(s,-1)+diag(s,1));
[X,I]=sort(diag(X)); x=(X+1)/2; wx=ab(1,2)*V(1,I)'.^2/4;
N=N-1; N1=N+1; N2=N+2;  y=cos((2*(N:-1:0)'+1)*pi/(2*N+2));
L=zeros(N1,N2);  y0=2;  iter=0;
while max(abs(y-y0))>eps
    L(:,1)=1;    L(:,2)=y;
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
    y0=y;    y=y0-L(:,N2)./Lp;  iter=iter+1;
end
cd=[ 1, 0, 0; -1, 0, 1; 0, 1,-1]*v;
t=(1+y)/2;  Wx=abs(det(cd(2:3,:)))*wx;  Wy=1./((1-y.^2).*Lp.^2)*(N2/N1)^2;
[tt,xx]=meshgrid(t,x); yy=tt.*xx;
X=cd(1,1)+cd(2,1)*xx+cd(3,1)*yy;    Y=cd(1,2)+cd(2,2)*xx+cd(3,2)*yy;
end
function area = findArea(p1,p2,p3)
u = p2 - p1;
v = p3 - p1;
biV = cross(u,v);
n = biV/norm(biV);
area = norm(biV)*0.5;
end
function out = indexNM(n,m)
mtot = 2*(n) + 1;
nvirtual = n + 1;
mvirtual = m + (mtot - 1)/2  + 1;
%Convert the n,m index to a position at a row vector Mnm
a1 = 1;
an = 2*(n-1)+1;
totN = (a1 + an)*(n)/2;
out = totN + mvirtual;
end
function out = Sb(n,y,m)
r = norm(y);
ir2 = 1/dot(y,y);
if n == 1
    out = 1/r;
    return;
end
%% Initialize Values
SValues(n) = SphericalHarmonics();
for i = 1:n
    SValues(i) = SphericalHarmonics(i);
end
%% Obtain N,N Values
SValues(1).addNNValue(1/r);
SValues(1).merge();
for i = 2:n
    Snn = SValues(i-1).MValuesPos(end);
    value = -Snn*(2*(i-2)+1)* (y(1) + complex(0,1)*y(2)) *ir2;
    SValues(i).addNNValue(value);
end
%% Obtain N,M Values
S10 = ir2*y(3)*SValues(1).MValuesPos(1);
SValues(2).addMValues(S10,1);
SValues(2).merge();
for i = 3:n
    obj = SValues(i);
    objAn2 = SValues(i-2);
    objAn1 = SValues(i-1);
    mEnd = 1+(obj.m-1)/2;
    for j = 1:mEnd-2
        MValue = ir2*(-((i-2) + (j-1))*((i-2) - (j-1))*objAn2.MValuesPos(j) + ...
            (2*(i-2) + 1)*y(3)*objAn1.MValuesPos(j) );
        obj.addMValues(MValue,j);
    end
    j = j + 1;
    MValue = ir2*(2*(i-2) + 1)*y(3)*objAn1.MValuesPos(j);
    obj.addMValues(MValue,j);
    obj.merge();
end
out = [SValues(:).values]';
if ~isempty(m)
    level = SValues(end).values';
    out = level(m);
end
end
