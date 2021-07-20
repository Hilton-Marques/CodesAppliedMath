function A = GAZ(B)
% GAZ(B): Set to zeros any element of the multivector B less than 1e-15
%
%See also gable.

% GABLE, Copyright (c) 1999, University of Amsterdam
% Copying, use and development for non-commercial purposes permitted.
%          All rights for commercial use reserved; for more information
%          contact Leo Dorst (leo@wins.uva.nl).
%
%          This software is unsupported.
A=B;
for i=1:8
  if abs(A.m(i)) < 1e-15
     A.m(i) = 0;
  end
end
