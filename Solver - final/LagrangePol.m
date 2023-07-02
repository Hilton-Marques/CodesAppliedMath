clear all;
no = 6;
x = linspace(0,1,no);
for f=1:no
pt1 = @(ksi) 1;
pt2 = @(ksi) 1;
for j = 1:f-1
pt1 = @(ksi) (ksi - x(j))/(x(f) - x(j)).*pt1(ksi)
end
for j = f+1:no
pt2 = @(ksi) (ksi - x(j))/(x(f) - x(j)).*pt2(ksi)
end
L = @(ksi) pt1(ksi)*pt2(ksi)
shape{f} = L;
end
shape{1}(0.1231);
shape{4}(0.6241);
shape{3}(0.6241);
indexat = @(expr,index) expr(index);

ksi = [2,3]
%f = diff(shape{1}([ksi ksi+0.00001] ))/0.00001
xksi = @(ksi) 3.*ksi;
yksi = @(ksi) 2.*ksi;
h = 0.00001;
r = @(ksi) [ xksi(ksi) ; yksi(ksi) ]
df = @(ksi) ((r(ksi+h)-r(ksi))) ./h;

n = @(ksi) [indexat(df(ksi),2);-indexat(df(ksi),1)] ./vecnorm(df(ksi));


