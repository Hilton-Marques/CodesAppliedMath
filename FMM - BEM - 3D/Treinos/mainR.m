clc;
clear;
close all;
%% Input Data
y = [0.0417, -0.0417, -0.1250];
n = 2;
R(n,y,1)
function out = R(n,y,flagAll)
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
    m = obj.m;
    mEnd = 1+(m-1)/2;
    for j = 1:mEnd-2
        MValue = (1/(((i-2) + (j-1) + 1)*(i-2 + 1 - (j-1))))*...
            ( (2*(i-2)+1)*y(3)*objAn1.MValuesPos(j) - ...
            r2*objAn2.MValuesPos(j) );
        obj.addMValues(MValue,j);
    end
    j = j +1;
    MValue = (1/(((i-2) + (j-1) + 1)*(i-2 + 1 - (j-1))))*...
        (2*(i-2)+1)*y(3)*objAn1.MValuesPos(j);
    obj.addMValues(MValue,j);
    obj.merge();
end
out = RValues(end).values;
if (flagAll == 1)
    out = [RValues(:).values];
end
end

