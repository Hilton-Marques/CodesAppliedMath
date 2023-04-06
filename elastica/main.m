clear
close(findall(0,'Type','figure'));
clc;

p0 = [0;0];

p1 = [1;3];
theta = pi/4.3
r = rot(-theta);
s = r*p1;
p2 = [2;-1];
p3 = s + p2;
pts = [p0,p1,p3,p2];
L3 = p3 - p1;
quiver(0,0,L3(1),L3(2));
drawQuad(pts);
L2 = norm(p2);
L1 = norm(p1);
fi = acos(dot(p2,p1)/(L2*L1));
sqrt(2.*L1.^2+L2.^2+(-2).*L1.*(L2.*cos(fi)+(-1).*L2.*cos(fi+(-1).* ...
  theta)+L1.*cos(theta)))
sqrt((L2+L1.*((-1).*cos(fi)+cos(fi+(-1).*theta))).^2+L1.^2.*(sin(fi)+( ...
  -1).*sin(fi+(-1).*theta)).^2)
norm(L3)
(dot(L3,p2)/(norm(L3)*norm(p2)))
(L2+(-1).*L1.*cos(fi)+L1.*cos(fi+(-1).*theta))/norm(L3)
keyboard
function drawQuad(pts)
   for i = 1:4
       pi = pts(:,i);
       pj = pts(:, mod(i,4)+1);
       line([pi(1),pj(1)],[pi(2),pj(2)],'linewidth',1);
   end
end



function r = rot(theta)
r = [[cos(theta),-sin(theta)];[sin(theta),cos(theta)]];
end