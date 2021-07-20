clc;
clear;
close all;
%% Input Data
%points = xlsread('Input','Planilha1');
tic
points = makeRand(100000);
N = size(points,1);
[~,I] = sort(points(:,1));
points = points(I,:);
convexHull = startConvexHull(points);
toc
% %% Plot
% hold on
% for i = 1:N
%     p = [points(i,1),points(i,2)];
%     plot(p(1),p(2),'o','Color','Red');
%     text(p(1),p(2),num2str(i),'HorizontalAlignment','Left')
% end
% convexHull.plot([0,0,0])
function convexHull = startConvexHull(points)
N = size(points,1);
if N <= 3
    convexHull = ConvexHull(points);
    %convexHull.plot([0,0,0])
    return
end
left_half = startConvexHull(points(1:fix(N/2),:));
right_half = startConvexHull(points(fix(N/2) + 1:N,:));
convexHull = merge(left_half,right_half);
%convexHull.plot([rand,rand,rand])
end
function convexHull = merge(left_half,right_half)
n1 = left_half.N;
n2 = right_half.N;
[~,iL] = max(left_half.points(:,1));
[~,iR] = min(right_half.points(:,1));
pLeft = left_half.points;
pRight = right_half.points;
%% Find upper Tangent
iLU = iL;
iRU = iR;
done = true;
while (done == true)
    done =  false;
    newIRU = mod(iRU,n2)+1;
    while(orient2D(pLeft(iLU,:),pRight(iRU,:),pRight(newIRU,:)) >= 0 )
        iRU = newIRU;
        newIRU = mod(iRU,n2)+1;
    end
    newILU = mod(n1 + iLU - 2,n1) + 1;
    while (orient2D(pRight(iRU,:),pLeft(iLU,:),pLeft(newILU ,:)) <= 0)
        iLU = newILU;
        newILU = mod(n1 + iLU - 2,n1)+1;
        done = true;
    end
end
%% Find low Tangent
iLL = iL;
iRL = iR;
done = true;
while (done == true)
    done =  false;
    newIRL = mod(n2 + iRL - 2,n2) +1 ;
    while(orient2D(pLeft(iLL,:),pRight(iRL,:),pRight(newIRL,:)) <= 0 )
        iRL = newIRL;
        newIRL = mod(n2 + iRL - 2,n2) + 1;
    end
    newILL = mod(iLL ,n1) + 1;
    while (orient2D(pRight(iRL,:),pLeft(iLL,:),pLeft(newILL ,:)) >= 0)
        iLL = newILL;
        newILL = mod(iLL ,n1)+1;
        done = true;
    end
end
%% Get index for new convexHUll clockwise
indexUpperLeft = 1:iLU;
indexRight = iRU:iRL;
if (iRU > iRL)
    indexRight = [iRU:n2,1:iRL];
end
indexLowLeft = iLL:n1;
if iLL == 1
    indexLowLeft = [];
end
convexHull = ConvexHull([pLeft(indexUpperLeft,:);...
    pRight(indexRight,:);pLeft(indexLowLeft,:)]);
end
function out = orient2D(a,b,c)
u = b - a;
v = c - b;
uPerp = [-u(2),u(1)];
out = dot(uPerp,v);
end
function points = makeRand(n)
points = zeros(n,2);
for i = 1:n
    points(i,:) = [rand,rand]; 
end
end