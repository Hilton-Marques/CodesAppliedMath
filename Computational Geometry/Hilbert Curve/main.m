clc;
clear;
close all;

order = 1;
n = 2^order;
total = n^2;

image = ones(300,400,3);
image(1:300,1:100,1) = 1;
imshow(image)


function out = hilbert(i)
 coord = [0,0;...
          0,1;...
          1,1;...
          1,0];
 out = coord(i,:);
 
end
