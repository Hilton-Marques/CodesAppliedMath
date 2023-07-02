% Clear workspace
clear
close(findall(0,'Type','figure'));
clc
%Set Initial State
n = 10;
initialState = makeTestRand(n);
%Number of iterations and threshold for Berner-Hut algorithm
hold on
pp = [initialState.position];
Points = [ [pp(1:2:end)]' [pp(2:2:end)]'];
draw_obj(Points);
tic 
Field(initialState,1);
toc

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




