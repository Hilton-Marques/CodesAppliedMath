clear
clc;
close all

Discrimant3(1,2,3)
Discrimant3(2,1,3)
Discrimant3(1,3,2)
Discrimant3(3,2,1)
Discrimant3(2,3,1)
Discrimant3(3,1,2)




function r = Discrimant3(a,b,c)
r = ((a-b)^(2/3) * (a - c)^(2/3) * (b - c)^(2/3));
end

function r = Discrimant2(a,b)
r = ((a-b)^2 )^1/2;
end