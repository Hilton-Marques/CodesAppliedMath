% Compute the error of a pose-pose constraint
% x1 3x1 vector (x,y,theta) of the first robot pose
% x2 3x1 vector (x,y,theta) of the second robot pose
% z 3x1 vector (x,y,theta) of the measurement
%
% You may use the functions v2t() and t2v() to compute
% a Homogeneous matrix out of a (x, y, theta) vector
% for computing the error.
%
% Output
% e 3x1 error of the constraint
% A 3x3 Jacobian wrt x1
% B 3x3 Jacobian wrt x2
function [e, A, B] = linearize_pose_pose_constraint(x1, x2, z)

% TODO compute the error and the Jacobians of the error

% Compute the homogenous transformations
X1 = v2t(x1);
X2 = v2t(x2);
Z = v2t(z);

% Compute the error for this constraint
%e = t2v(inv(Z)*(inv(X1)*X2));
E = Z\(X1\X2);
 global flag;

if flag
    e = log(E);
    A = getJacobianInv(e)*(-inv(getAdj(X2))*(getAdj(X1)));
    B = getJacobianInv(e);
    return;
end
e = t2v(E);


% Note: The code below only requires the computations of relative
% measurement and pose relative angles. However, this will make the
% code look a little involved. Thus, in order to compute the Jacobians blocks
% each of the necesarry components for the computation are set separately.

% Set the measurment rotation matrix
Rij = Z(1:3, 1:3);

% Set pose X1 rotation matrix
Ri = X1(1:3, 1:3);

% Compute the pose X1 relative angle
theta_i = atan2(Ri(2, 1),Ri(1, 1));

% Compute the measurement relative angle
theta_ij = atan2(Rij(2, 1),Rij(1, 1));

% Set pose X1 coordinates
xi = x1(1);
yi = x1(2);

% Set pose X2 coordinates
xj = x2(1);
yj = x2(2);

% Compute eij partial derivative with respect to pose x1
eij_xi = [-cos(theta_i)*cos(theta_ij)+sin(theta_i)*sin(theta_ij);
    cos(theta_i)*sin(theta_ij)+sin(theta_i)*cos(theta_ij);
    0
    ];


eij_yi = [-sin(theta_i)*cos(theta_ij)-cos(theta_i)*sin(theta_ij);
    sin(theta_i)*sin(theta_ij)-cos(theta_i)*cos(theta_ij);
    0
    ];

eij_theta_i = [-(xj - xi)*(sin(theta_i)*cos(theta_ij)+cos(theta_i)*sin(theta_ij))+(yj - yi)*(cos(theta_i)*cos(theta_ij)-sin(theta_i)*sin(theta_ij));
    (xj - xi)*(sin(theta_i)*sin(theta_ij) - cos(theta_i)*cos(theta_ij))-(yj - yi)*(cos(theta_i)*sin(theta_ij)+sin(theta_i)*cos(theta_ij));
    -1
    ];

% Construct Jacobian block A with respect to x1
A = [eij_xi, eij_yi, eij_theta_i];

% Compute the partial derivative with respect to pose x2
eij_xj = [cos(theta_i)*cos(theta_ij)-sin(theta_i)*sin(theta_ij);
    -cos(theta_i)*sin(theta_ij)-sin(theta_i)*cos(theta_ij);
    0
    ];
eij_yj = [sin(theta_i)*cos(theta_ij)+cos(theta_i)*sin(theta_ij);
    -sin(theta_i)*sin(theta_ij)+cos(theta_i)*cos(theta_ij);
    0
    ];

eij_theta_j = [0; 0; 1];

% Construct Jacobian block B with respect to x2

B = [eij_xj, eij_yj, eij_theta_j];
end
function Jinv = getJacobianInv(v)
theta = v(3);
v1 = v(1);
v2 = v(2);
cot_tan_2 = cos(theta/2)/sin(theta/2);
Jinv = [(1/2).*theta.*cot_tan_2,(-1/2).*theta,(1/2).*(v2+v1.*(2.*theta.^(-1)+(-1) ...
    .*cot_tan_2));(1/2).*theta,(1/2).*theta.*cot_tan_2,(-1/2).*v1+v2.* ...
    theta.^(-1)+(-1/2).*v2.*cot_tan_2;0,0,1];
end
function adj = getAdj(P)
adj = [P(1:2,1:2),-[-P(2,3);P(1,3)];[0,0,1]];
end
function Q = exp(v,P)
theta = v(3);
t = v(1:2);

%Rotation part
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

%Translation part
V = eye(2);
if (theta ~= 0)
    V = [sin(theta), -(1-cos(theta)); 1-cos(theta), sin(theta)] .* (1/theta);
end
t = V * t;

Q = P*[R, t; 0 0 1];
end
function v = log(C)
w = atan2(C(2,1), C(1,1));
w_sqr = w*w;
if (w_sqr <  1e-14)
    %taylor approximation
    a = 1 - (1/6)*w_sqr;
    b = 0.5 * w - w*w_sqr*(1/24);
else
    a = sin(w)/w;
    b = (1-cos(w))/w;
end
invV = [a, b; -b, a] .* 1./(a^2+b^2);
t = invV * C(1:2,end);
v = [t;w];
end
