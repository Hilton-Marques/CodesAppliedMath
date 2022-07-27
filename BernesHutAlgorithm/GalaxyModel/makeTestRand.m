function initialState = makeTestRand(starNumber)
% This function randomly creates a cluster of stars.
%   intialState is an array of structs whose members are mass, position,
%   velocity, acceleration.
%   rand:  Uniformly distributed pseudorandom numbers.
%   randn: Normally distributed pseudorandom numbers.

axisLim = 5;
initialState(starNumber).mass = 0;
for s = 1:starNumber
    initialState(s).mass = 1;
    initialState(s).position = 0.5* axisLim * [rand,rand,rand];
    initialState(s).velocity = 0*[rand,rand,rand];
    initialState(s).acceleration = [0,0,0];
end

end

