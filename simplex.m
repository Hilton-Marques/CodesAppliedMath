clc;
clear;
close all;

A = [1,2,3,0;-1,2,6,0;0,4,9,0;0,0,3,1];

I = eye(4,3);
newA = [A,I];
c = [zeros(4,1);ones(3,1)];
basis = 8-3 : 8;
basis = [5,6,7,4];
b = [3;2;5;1];
n = 5;
newA = [1,1,1,0,0;1,0,0,1,0;0,1,0,0,1];
b = [20;16;10];
c = [-2;-3;0;0;0];
basis = [3,4,5];
x = solveLP(newA,b,c,basis)
[xOther,fval] = linprog(c,-eye(n,n),zeros(n,1),newA,b);

function x = solveLP(A,b,c,basis)
At = A';
x = zeros(size(c,1),1);
d = x;
dAux = x;
B = A(:,basis);
x(basis) = b; % xpos at starting corner
cost = c(basis)'*x(basis);     % cost at starting corner
n = size(A,2);
for iter = 1:100
   % d = zeros(size(c,1),1);
    %dAux = zeros(size(c,1),1);
    y = transpose(inv(B))*c(basis);           % this y may not be feasible
    %nonBasics = setdiff(1:n,basis);
    rmin = zeros(n,1);
    for in = 1:n
        rmin_i = c(in) - At(in,:)*y;
        rmin(in) = rmin_i;
%         if rmin < 0
%             break
%         end
    end
    [rmin,in] = min(rmin);
    %[rmin,in] = min(c - A'*y); % minimum r and its index in
    if rmin >= -.00000001      % optimality is reached, r>=0
        check = A*x - b;
        break;                 % current x and y are optimal
    end
    
    d(basis) = B\A(:,in);  % decrease in x from 1 unit of xin
    dAux(basis) =  -d(basis);
    dAux(in) = 1;
    A*dAux
    xb = x(basis);
    db = d(basis);
    ids = find(db > 0 );
    tetas = xb(ids)./db(ids);
    [teta,~] = min(tetas);
    %Blands rule
    if teta == 0
        outs = find(teta==xb);
    else
        outs = ids(find(teta==tetas));
    end
    out = outs(1);
    for i = 2:size(outs)
        outTemp = outs(i);
        if basis(outTemp) < basis(out)
            out = outTemp;
        end       
    end
    if d(out) == .000001  % out = index of first x to reach 0
        break;      % break when that edge is extremely short
    end
    cost = cost + teta*rmin;  % lower cost at end of step
    x(basis) = x(basis) - teta*d(basis);   % update old x
    x(in) = teta;      % find new positive component of x
    check = A*x - b
    basis(out) = in;      % replace old index by new in basis
    basis = sort(basis);
    B = A(:,basis);
end

end