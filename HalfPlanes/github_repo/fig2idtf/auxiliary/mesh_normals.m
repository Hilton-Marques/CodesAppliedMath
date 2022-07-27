function [normals] = mesh_normals(points,faces)

% MESH_NORMALS   compute mesh normals
%
%   compute mesh normals using builtin matlab function
%
%   SYNTAX
%       [NORMALS] = MESH_NORMALS(POINTS,FACES)
%
%   Created by Alexandre Gramfort on 2007-11-27.
%   Copyright (c) 2007 Alexandre Gramfort. All rights reserved.

% $Id: mesh_normals.m 26 2008-11-03 10:14:12Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2008-11-03 11:14:12 +0100 (Lun, 03 nov 2008) $
% $Revision: 26 $

me = 'MESH_NORMALS';

if nargin == 0
    eval(['help ',lower(me)])
    return
end

hf = figure('Visible','off');
hp = patch('vertices',points,'faces',faces);
normals = get(hp,'Vertices');
close(hf);
%Make the normals unit norm
norms = sqrt(sum(normals.*conj(normals),2));
gidx = find(norms);
normals(gidx,:) = - normals(gidx,:) ./ repmat(norms(gidx),1,3);

% n = size(points,1);
% normals = zeros(n,3);
% for i=1:size(faces,1)
%    A=faces(i,:);
%    B=points(A',:);
%    a=B(1,:);
%    b=B(2,:);
%    c=B(3,:);
%    
%    ab=a-b;
%    ac=a-c;
%    
%    norm=cross(ab,ac);
%    normMag=sqrt(sum(norm.^2));
%    normals(A,:)= normals(A,:) + ones(3,1) .* norm/normMag;
% end
% norms = sqrt(sum(normals.*conj(normals),2));
% gidx = find(norms);
% normals(gidx,:) = normals(gidx,:) ./ repmat(norms(gidx),1,3);% normals = repmat(normals,3,1);