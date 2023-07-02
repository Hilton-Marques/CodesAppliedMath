clc;
clear;
close all;
%% Input Data
n = 3;
y = [1.2,2.2,3.5];
Sb(n,y)
function out = Sb(n,y)
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
    m = obj.m;
    mEnd = 1+(m-1)/2;
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
out = SValues(end).values';
end