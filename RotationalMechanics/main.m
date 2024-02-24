% Motion of a book (rectangular rigid body) tossed in the air

% Written by Alyssa Novelia (4/24/2013)

% Revised by Daniel Kawano (12/21/2017)

clear all;
close all;
clc;

%--------------------------------------------------------------------------

% Physical and simulation parameters:

% Set up a rectangular object with length L, width w, and thickness b.
% Length is associated with the e1 corotational basis vector, width with 
% e2, and thickness with e3.

L = 25/100;     % m
w = 18/100;     % m
b = 2.5/100;    % m
m = 0.7;        % kg
g = 0;        % m/s^2

% Calculate the principal moments of inertia:

lambda1 = (m/12)*(w^2 + b^2); 	% kg-m^2
lambda2 = (m/12)*(L^2 + b^2); 	% kg-m^2
lambda3 = (m/12)*(w^2 + L^2);  	% kg-m^2

% Initial mass center velocity:

xdot0 = [0, 0, 0];          % [x1dot(0), x2dot(0), x3dot(0)], m/s

% Initial mass center position:

x0 = [0, 0, 0];             % [x1(0), x2(0), x3(0)], m

% Initial angular velocity: 

omega0 = [0.1, 20, 0.1];   	% [omega1(0), omega2(0), omega3(0)], rad/s 

% Initial orientation/Euler angles:

Q0 = [0, 0, 0];             % [psi(0), theta(0), phi(0)], rad

% Simulation parameters:

dt = 0.005;
tf = 1;
tsim = [0 : dt : tf]';

tol = 1e-8;
options = odeset('abstol', tol, 'reltol', tol);

%--------------------------------------------------------------------------

% Integrate the equations of motion:

[T, Y] = ode45(@EOM, tsim, [xdot0, x0, omega0, Q0], options, m, g, lambda1, lambda2, lambda3);

% Extract the position and orientation solutions:

x1 = Y(:,4);        % m
x2 = Y(:,5);        % m
x3 = Y(:,6);        % m
psi = Y(:,10);      % rad
theta = Y(:,11);    % rad
phi = Y(:,12);      % rad

%--------------------------------------------------------------------------

% Compute the rotation tensor components and the evolution of the
% corotational basis vectors:

% The rotation is parameterized using 3-2-1 Euler angles:

e1seq = zeros(3,length(phi));
e2seq = zeros(3,length(phi));
e3seq = zeros(3,length(phi));

for iii = 1:length(phi)
    
    PSI = psi(iii);
    THE = theta(iii);
    PHI = phi(iii);
    
    Qmat = [cos(PSI) -sin(PSI) 0; sin(PSI) cos(PSI) 0; 0 0 1]*...
        [cos(THE) 0 sin(THE); 0 1 0; -sin(THE) 0 cos(THE)]*...
        [1 0 0; 0 cos(PHI) -sin(PHI); 0 sin(PHI) cos(PHI)];
    
    e1seq(:,iii) = Qmat(:,1);
    e2seq(:,iii) = Qmat(:,2);
    e3seq(:,iii) = Qmat(:,3);

end

%--------------------------------------------------------------------------

% Plot the position of the mass center and the book's orientation over
% time:

% figure;
% set(gcf, 'color', 'w', 'name', 'Mass center position and orientation over time');
% 
% subplot(1,2,1);
% plot(T, x1, T, x2, T, x3, 'linewidth', 2);
% title('Mass center position over time');
% xlabel('Time (s)');
% ylabel('Mass center position (m)');
% legend('\itx\rm_1', '\itx\rm_2', '\itx\rm_3');
% 
% subplot(1,2,2);
% plot(T, psi*(180/pi), T, theta*(180/pi), T, phi*(180/pi), 'linewidth', 2);
% title('Orientation over time');
% xlabel('Time (s)');
% ylabel('Angle (deg)');
% legend('\psi', '\theta', '\phi');

%--------------------------------------------------------------------------

% Plot the evolution of the corotational basis by tracing the trajectories
% of the basis vectors' tips on the unit sphere:

% figure;
% set(gcf, 'color', 'w', 'name', 'Evolution of the corotational basis');
% 
% [X, Y, Z] = sphere;
% sphere_obj = surf(X, Y, Z);
% set(sphere_obj, 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1);
% 
% title({'Trajectories traced by the corotational basis on the unit sphere:'; 'e_1 = blue, e_2 = red, e_3 = black'});
% xlabel('\bfE\rm_1'); 
% ylabel('\bfE\rm_2'); 
% zlabel('\bfE\rm_3  ', 'rotation', 0);
% axis equal;
% hold on;
% 
% % Draw the basis vectors:
% 
% e1v = quiver3(0,0,0, e1seq(1,1), e1seq(2,1), e1seq(3,1), 'b', 'LineWidth', 2, 'AutoScale', 'off');
% e2v = quiver3(0,0,0, e2seq(1,1), e2seq(2,1), e2seq(3,1), 'r', 'LineWidth', 2, 'AutoScale', 'off');
% e3v = quiver3(0,0,0, e3seq(1,1), e3seq(2,1), e3seq(3,1), 'k', 'LineWidth', 2, 'AutoScale', 'off');
% 
% % Draw their tips' trajectories:
% 
% e1p = line(e1seq(1,:), e1seq(2,:), e1seq(3,:), 'Color', 'b', 'LineWidth', 1.5);
% e2p = line(e2seq(1,:), e2seq(2,:), e2seq(3,:), 'Color', 'r', 'LineWidth', 1.5);
% e3p = line(e3seq(1,:), e3seq(2,:), e3seq(3,:), 'Color', 'k', 'LineWidth', 1.5);

%--------------------------------------------------------------------------

% Animating the motion of the tossed book:

% Set up the animation window:

figure;
set(gcf, 'color', 'w', 'name', 'Animation of the tossed book');
axis equal;
xlim([min(x1)-L, max(x1)+L]);
ylim([min(x2)-L, max(x2)+L]);
zlim([min(x3)-L, max(x3)+L]);
xlabel('\itx\rm_1 (m)'); 
ylabel('\itx\rm_2 (m)'); 
zlabel('\itx\rm_3 (m)         ', 'rotation', 0);
view([135 30]);
grid on;

% We need to animate 6 planes to form a rectangular prism. Track 8 material 
% points, chosen to be the vertices of the prism:

rCM = [x1, x2, x3]';

r1 = rCM + 0.5*L*e1seq + 0.5*w*e2seq + 0.5*b*e3seq;
r2 = rCM - 0.5*L*e1seq + 0.5*w*e2seq + 0.5*b*e3seq;
r3 = rCM - 0.5*L*e1seq - 0.5*w*e2seq + 0.5*b*e3seq;
r4 = rCM + 0.5*L*e1seq - 0.5*w*e2seq + 0.5*b*e3seq;
r5 = rCM + 0.5*L*e1seq + 0.5*w*e2seq - 0.5*b*e3seq;
r6 = rCM - 0.5*L*e1seq + 0.5*w*e2seq - 0.5*b*e3seq;
r7 = rCM - 0.5*L*e1seq - 0.5*w*e2seq - 0.5*b*e3seq;
r8 = rCM + 0.5*L*e1seq - 0.5*w*e2seq - 0.5*b*e3seq;

% Calculate the vertices locations over time for the 6 planes:

vertices1_x = [r1(1,:); r2(1,:); r3(1,:); r4(1,:)];
vertices1_y = [r1(2,:); r2(2,:); r3(2,:); r4(2,:)];
vertices1_z = [r1(3,:); r2(3,:); r3(3,:); r4(3,:)];
surf_1 = patch(vertices1_x(:,1), vertices1_y(:,1), vertices1_z(:,1), 'FaceColor', 'g', 'FaceAlpha', 0.5);

vertices2_x = [r5(1,:); r6(1,:); r7(1,:); r8(1,:)];
vertices2_y = [r5(2,:); r6(2,:); r7(2,:); r8(2,:)];
vertices2_z = [r5(3,:); r6(3,:); r7(3,:); r8(3,:)];
surf_2 = patch(vertices2_x(:,1), vertices2_y(:,1), vertices2_z(:,1), 'FaceColor', 'g', 'FaceAlpha', 0.5);

vertices3_x = [r3(1,:); r4(1,:); r8(1,:); r7(1,:)];
vertices3_y = [r3(2,:); r4(2,:); r8(2,:); r7(2,:)];
vertices3_z = [r3(3,:); r4(3,:); r8(3,:); r7(3,:)];
surf_3 = patch(vertices3_x(:,1), vertices3_y(:,1), vertices3_z(:,1), 'FaceColor', 'g', 'FaceAlpha', 0.5);

vertices4_x = [r1(1,:); r2(1,:); r6(1,:); r5(1,:)];
vertices4_y = [r1(2,:); r2(2,:); r6(2,:); r5(2,:)];
vertices4_z = [r1(3,:); r2(3,:); r6(3,:); r5(3,:)];
surf_4 = patch(vertices4_x(:,1), vertices4_y(:,1), vertices4_z(:,1), 'FaceColor', 'g', 'FaceAlpha', 0.5);

vertices5_x = [r1(1,:); r4(1,:); r8(1,:); r5(1,:)];
vertices5_y = [r1(2,:); r4(2,:); r8(2,:); r5(2,:)];
vertices5_z = [r1(3,:); r4(3,:); r8(3,:); r5(3,:)];
surf_5 = patch(vertices5_x(:,1), vertices5_y(:,1), vertices5_z(:,1), 'FaceColor', 'g', 'FaceAlpha', 0.5);

vertices6_x = [r2(1,:); r3(1,:); r7(1,:); r6(1,:)];
vertices6_y = [r2(2,:); r3(2,:); r7(2,:); r6(2,:)];
vertices6_z = [r2(3,:); r3(3,:); r7(3,:); r6(3,:)];
surf_6 = patch(vertices6_x(:,1), vertices6_y(:,1), vertices6_z(:,1), 'FaceColor', 'g', 'FaceAlpha', 0.5);

% Highlight one of the vertices:

point6 = line(r6(1,1), r6(2,1), r6(3,1), 'marker', 'o', 'markerfacecolor', 'b');

% Animate the tossed book:

% animation = VideoWriter('tossed-book.avi');
% animation.FrameRate = 1/dt/4;       
% open(animation);
im = {};
for jjj = 1:length(e1seq)
    
    set(surf_1, 'xdata', vertices1_x(:,jjj), 'ydata', vertices1_y(:,jjj), 'zdata', vertices1_z(:,jjj));
    set(surf_2, 'xdata', vertices2_x(:,jjj), 'ydata', vertices2_y(:,jjj), 'zdata', vertices2_z(:,jjj));
    set(surf_3, 'xdata', vertices3_x(:,jjj), 'ydata', vertices3_y(:,jjj), 'zdata', vertices3_z(:,jjj));
    set(surf_4, 'xdata', vertices4_x(:,jjj), 'ydata', vertices4_y(:,jjj), 'zdata', vertices4_z(:,jjj));
    set(surf_5, 'xdata', vertices5_x(:,jjj), 'ydata', vertices5_y(:,jjj), 'zdata', vertices5_z(:,jjj));
    set(surf_6, 'xdata', vertices6_x(:,jjj), 'ydata', vertices6_y(:,jjj), 'zdata', vertices6_z(:,jjj));
    set(point6, 'xdata', r6(1,jjj), 'ydata', r6(2,jjj), 'zdata', r6(3,jjj));
    drawnow;
    im = get(im);
    %writeVideo(animation, getframe(gcf));
end
save(im, filename="cube.gif");
% close(animation);

%==========================================================================

%  Specify the equations of motion:

function dY = EOM(T, Y, m, g, lambda1, lambda2, lambda3)

dY = zeros(12,1);
    
x1dot = Y(1);
x2dot = Y(2);
x3dot = Y(3);
x1 = Y(4);
x2 = Y(5);
x3 = Y(6);
omega1 = Y(7);
omega2 = Y(8);
omega3 = Y(9);
psi = Y(10);
theta = Y(11);
phi = Y(12);

M = blkdiag(m*eye(3), ...
            eye(3), ...
            blkdiag(lambda1, lambda2, lambda3), ...
            [-sin(theta), 0, 1;
             cos(theta)*sin(phi), cos(phi), 0;
             cos(theta)*cos(phi), -sin(phi), 0]);

f = [0;
    0;
    -m*g;
    x1dot;
    x2dot;
    x3dot;
    -(lambda3-lambda2)*omega2*omega3;
    -(lambda1-lambda3)*omega1*omega3;
    -(lambda2-lambda1)*omega1*omega2;
    omega1;
    omega2;
    omega3];

dY = M\f;
end

 function save(im,options)
            arguments  
                im
                options.filename = filename
                options.slow = false
                options.repeat = false
            end
            
            for idx = 1:length(im)
                [A,map] = rgb2ind(im{idx},256);
                if idx == 1
                    imwrite(A,map, options.filename,'gif','LoopCount',Inf,'DelayTime',1.5);
                elseif (idx > 1) && (idx < length(im))
                    if (options.slow)
                        imwrite(A,map, options.filename,'gif','WriteMode','append','DelayTime',0.16);
                    else
                        imwrite(A,map, options.filename,'gif','WriteMode','append','DelayTime',0.0278);
                    end
                else
                    if ~(options.repeat)
                        imwrite(A,map, options.filename,'gif','WriteMode','append','DelayTime',655);
                    else
                        imwrite(A,map, options.filename,'gif','WriteMode','append','DelayTime', 1);
                    end
                end
            end
 end

 function im = get(im)
 %frame = getframe(this.figure);
 frame = getframe(gcf);
 im{end+1} = frame2im(frame);
 end
