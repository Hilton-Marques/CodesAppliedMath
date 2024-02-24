%{
%Copyright (c) 2023 Hilton-Marques <https://my.github.com/Hilton-Marques>
%
%Created Date: Wednesday, September 27th 2023, 10:37:57 pm
%Author: Hilton-Marques
%
%Description: A program to study affine transformations
%HISTORY:
%Date      	By	Comments
%----------	---	----------------------------------------------------------
%}
clc;
clear;
close all;
addpath(genpath("../../my_libs/"));

Circle2Ellipse();



%----------------- Circle2Ellipse -----------------%
function Circle2Ellipse()
 theta = linspace(0,2*pi, 100);
 O = [cos(theta); sin(theta)];
 obj = AffineGeometry("circle2ellipse.gif");
 T = [[1,sqrt(2)/2];[sqrt(2)/2,2]];
 obj.TransformPolygon(T, O);
 obj.save(repeat=true);
 obj.exportFrame();
end

%-----------------Affine Triangle-----------------%
function AffineTriangle()
obj = AffineGeometry("affine-triangle.gif");
r = [1;1];
q = [-2;1];
p = [4;4];
r = p + r;
q = r + q;
obj.TransformTriangle(p,q,r);
obj.save();
obj.exportFrame();
end