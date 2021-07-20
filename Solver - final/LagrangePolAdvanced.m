clc;
clear all;
syms ksi
no = 5;
x = linspace(0,1,no);
f1 = 1;
for j = 1:no-1 
for i = 1: no-1
    if j == f1
        M{i,f1}= @(ksi) ksi^i-1;
    else
        M{i,j} = x(j).^i-1
    end
end
end
detX = @(ksi) cellfun(@det,M)
detX(0.3)

        
