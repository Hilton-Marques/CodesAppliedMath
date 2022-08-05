% Clear workspace
clear
%close(findall(0,'Type','figure'));
clc


p = gcp();
tic
n = 2;
% To request multiple evaluations, use a loop.
for idx = 1:n
  f(idx) = parfeval(p,@magic,1,idx); % Square size determined by idx
end
% Collect the results as they become available.
magicResults = cell(1,10);
for idx = 1:n
  % fetchNext blocks until next results are available.
  fprintf('hello')
  [completedIdx,value] = fetchNext(f);
  magicResults{completedIdx} = value;
  fprintf('Got result with index: %d.\n', completedIdx);
end
toc