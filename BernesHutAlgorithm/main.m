% Clear workspace
clear
close(findall(0,'Type','figure'));
clc
%Set Initial State
n = 500;
initialState = makeTestRand(n);
%Number of iterations and threshold for Berner-Hut algorithm
max_iter = 50;
dt = 0.01;
theta = 0.5;
for i = 1:max_iter
    %find convex hull
    pp = [initialState.position];
    Points = [ [pp(1:2:end)]' [pp(2:2:end)]'];
    boundary = [min(Points,[],1) max(Points,[],1)];
    cellCentroid = mean([boundary(1:2); boundary(3:4)], 1);
    s = max((boundary(3) - boundary(1)),(boundary(4) - boundary(2)));
    %Show bodies
    clf;
    hold on
    axis([cellCentroid(1,1) - s/2,cellCentroid(1,1) + s/2 ,...
        cellCentroid(1,2) - s/2, cellCentroid(1,2) + s/2]);
    draw_obj(Points);
    pause(0.0001)
    %Build tree
    capacity = 5; %max capacity for each leaf
    root = Leaf(s,cellCentroid,capacity);
    for j = 1:length(initialState)
        points = initialState(j).position;
        root.insertBody(initialState(j));
    end
    %Update bodies
    verlet(initialState,root,dt,theta);
end
function verlet(initialState,root,dt,theta)
for i = 1: length(initialState)
    initialState(i).accelerationOld = initialState(i).acceleration;
    updateAccelerationOn(initialState(i),root,theta);
    initialState(i).advance(dt);
end
end
function updateAccelerationOn(body,root,theta)
% %Plot what is happening
% [h1,h2] = root.showCurrent(body);
% delete([h1,h2]);
rSqr = root.distSqr(body);
if (~root.isDivided)
    root.updateAccelerationOn(body,rSqr);
elseif ((root.s)^2/rSqr < (theta^2) )
    root.updateAccelerationOn(body,rSqr);
else
    for i = 1:4
        child = root.children{i,1};
        if (child.nBodies > 0)
            updateAccelerationOn(body,child,theta);
        end
    end
end
end


function initialState = makeTestRand(starNumber)
% This function randomly creates a cluster of stars.
%   intialState is an array of structs whose members are mass, position,
%   velocity, acceleration.
%   rand:  Uniformly distributed pseudorandom numbers.
%   randn: Normally distributed pseudorandom numbers.

axisLim = 5;
initialStateTemp(starNumber).mass = 0;
initialState(starNumber) = Body();
for s = 1:starNumber
    if (s==1)
        initialStateTemp(s).mass = 1;
        initialStateTemp(s).position = 0.5* axisLim * [rand,rand];
        initialStateTemp(s).velocity = 0*[rand,rand];
        initialStateTemp(s).acceleration = [0,0];
        initialState(s) = Body(initialStateTemp(s));
        
    else
        initialStateTemp(s).mass = 1;
        initialStateTemp(s).position = 0.5* axisLim * [rand,rand];
        initialStateTemp(s).velocity = 0*[rand,rand];
        initialStateTemp(s).acceleration = [0,0];
        initialState(s) = Body(initialStateTemp(s));
    end
end

end
function draw_obj(points)
plot(points(:,1),points(:,2),'o');
% for i = 1:size(points,1)
%     text(points(i,1),points(i,2),num2str(i),'HorizontalAlignment','left');
% end
end
