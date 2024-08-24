%{
%Copyright (c) 2023 Hilton-Marques <https://my.github.com/Hilton-Marques>
%
%Created Date: Thursday, November 23rd 2023, 10:42:28 am
%Author: Hilton-Marques
%
%Description: A  class laboratory to study rotational mechanics
%HISTORY:
%Date      	By	Comments
%----------	---	----------------------------------------------------------
%}

classdef IMU < Scene
    properties
        m_m % prism mass
        m_J % inertia tensor
        m_R % orientation
        m_scale %
        m_w %angular velocity
        m_dt
        m_h
        m_accel %acceleration
        m_mag %magnetometer
    end
    methods
        function this = IMU(l,w,b, wt, dt,accel,magno)
            this = this@Scene('IMU.gif','SE3');
            this.m_scale = diag([l,w,b]);
            this.m_R = SO3;
            this.m_R.setToIdentity();
            this.m_w = wt;
            this.m_dt = dt;
            this.m_accel = accel;
            this.m_mag = magno;
            this.BuildInertialTensor(l,w,b);
            %this.Init();
            axis("tight");
        end

        function Init(this)
            view(135,30);
            T = this.m_R.m_data';
            T(4,4) = 1.0;
            [h, pts] = this.drawRobot3D(T,scale=this.m_scale);
            this.setMargin(0.15);
            this.setBB(pts');
            this.get();
            this.m_h = h;
        end

        function BuildInertialTensor(this, l, w, b)
            lambda1 = (this.m_m/12) * (w^2 + b^2);  % kg-m^2
            lambda2 = (this.m_m/12) * (l^2 + b^2);  % kg-m^2
            lambda3 = (this.m_m/12) * (w^2 + l^2);  % kg-m^2
            this.m_J = diag([lambda1, lambda2, lambda3]);
        end

        function rots = Solver(this, plot, exact)
            if nargin == 1
                plot = false;
                exact = [];
            end
            if nargin == 2
                exact = [];
            end
            this.m_R.setToIdentity();
            n = size(this.m_w,1);
            dt = this.m_dt;
            rots = cell(1,n);
            rots{1} = this.m_R;
            quat_orientation = quaternion([1 0 0 0]);            
            R = eye(3);
            for t = 1:n-1
                wd = this.m_w(t,:);
                this.m_R = this.m_R + (wd * dt); % update
%                 R = R*expm(vec2skew(wd * dt));
%                 quat_orientation = quat_orientation * quaternion(wd*dt,"rotvec");
%                 quat_orientation_rotm = quat2rotm(quat_orientation);
                rots{t+1} = this.m_R;
                if (plot)
                    T = this.m_R.m_data;
                    T(4,4) = 1.0;
                    delete(this.m_h);
                    [this.m_h,pts] = this.drawRobot3D(T,scale=this.m_scale);
                    if (~isempty(exact))         
                        exact_rot = quat2rotm(exact(t+1));                        
                        exact_rot(4,4) = 1.0;                        
                        [h,pts] = this.drawRobot3D(exact_rot,...
                                                   scale=this.m_scale,...
                                                   color=this.m_blue);
                        this.m_h = [this.m_h, h];
                    end
                    %show angular velocity
                    scale = 0.3;
                    w_frame = scale*this.m_R.m_data*wd';                    
                    h = this.arrow([0;0;0], w_frame,color=this.m_red);
                    %pause(0.01);
                    %drawnow;
                    this.get();
                    delete(h);
                end
            end
            this.save(filename="raw_imu.gif");
        end

        function rots = SolverCMF(this, plot, exact)
            if nargin == 1
                plot = false;
                exact = [];
            end
            if nargin == 2
                exact = [];
            end
            this.m_R.setToIdentity();
            n = size(this.m_w,1);
            dt = this.m_dt;
            rots = cell(1,n);
            rots{1} = this.m_R;
            quat_orientation = quaternion([1 0 0 0]);            
            R = eye(3);
            bias = zeros(n,3);
            ki = 0.2; kp = 1.0;
            g0 = exact.g0;
            b0 = exact.B0;
            for t = 1:n-1
                wd = this.m_w(t,:);
                inv_R = this.m_R.m_data';
                R = this.m_R.m_data;
                sigma_r = cross(this.m_accel(t,:), g0*R) + ...
                          cross(this.m_mag(t,:), b0*R);                
                wd = wd - bias(t,:) + kp * sigma_r; 
                this.m_R = this.m_R + (wd * dt);
                bias(t+1,:) = bias(t,:) - ki * sigma_r * dt;
%                 R = R*expm(vec2skew(wd * dt));
%                 quat_orientation = quat_orientation * quaternion(wd*dt,"rotvec");
%                 quat_orientation_rotm = quat2rotm(quat_orientation);
                rots{t+1} = this.m_R;
                if (plot)
                    T = this.m_R.m_data;
                    T(4,4) = 1.0;
                    delete(this.m_h);
                    [this.m_h,pts] = this.drawRobot3D(T,scale=this.m_scale);
                    if (~isempty(exact))         
                        exact_rot = quat2rotm(exact.orientation(t+1));                        
                        exact_rot(4,4) = 1.0;                        
                        [h,pts] = this.drawRobot3D(exact_rot,...
                                                   scale=this.m_scale,...
                                                   color=this.m_blue);
                        this.m_h = [this.m_h, h];
                    end
                    %show angular velocity
                    scale = 0.3;
                    w_or = this.m_w(t,:);
                    w_frame_or = scale*this.m_R.m_data*w_or';        
                    h_or = this.arrow([0;0;0], w_frame_or,color=this.m_red);
                    
                    w_frame_curr = scale*this.m_R.m_data*wd';        
                    h_curr = this.arrow([0;0;0], w_frame_curr,color=this.m_blue);
                    %pause(0.01);
                    %drawnow;
                    this.get();
                    delete(h_or);
                    delete(h_curr);
                end
            end
            this.save(filename="imu_cmf.gif");
        end

        function ShowBias(this)
            filename = "white_noise_with_envelope";
            noise_white = true;
            noise_bias = false;
            only_bias = false;
            with_envelope = true;
            n = 10;           
            std = 0.05;
            x = linspace(0, 1, n);
            dt = x(end)/n;
            true_value = 9.81;
            y = 0*x + true_value;
            
            %add noise
            if (noise_white)
            y = y + normrnd(0, std, 1,n);
            end

            %add bias
            bias_std = 0.1*sqrt(dt);
            bias = 0*x;
            for i = 1:n-1
                bias(i+1) = bias(i) + normrnd(0, bias_std);
            end
            if noise_bias
                y = y + bias;
            end
            if only_bias
                y = bias;
            end

                        %plots
            fig = figure;
            ax = axes('Parent', fig); % Create axes in the figure
            set(gcf,'color','white');
            xaxis(1);
            
            xticks(0:0.1:1.0);
            if (~noise_bias)
                yaxis([9.5, 10.0])
                yticks([9.50, 9.81, 10]); % Specify the ticks you want on the y-axis
                yticklabels({'9.50', '9.81', '10.00'}); % Custom labels for each tick
            else
                yaxis([min(y) - 0.25, max(y) + 0.25]);
            end
            grid on
            set(ax,'XMinorGrid','off');
            set(ax,'YMinorGrid','off');
            set(ax,'FontSize',12);
            set(ax,'XMinorTick','off');
            set(ax,'YMinorTick','on');            
            xlabel(strcat('$', "t(s)", '$'),'interpreter','latex');            
            ylabel(strcat('$', "\tilde{a}_{B}^{z}(m/s^2)", '$'),'interpreter','latex');
            hold(ax, 'on'); % Hold on to the current plot            
            if ~with_envelope
                for i = 1:n
                    plot(ax, x(1:i), y(1:i),'color','black','linewidth',1);
                    this.get(fig);
                end
            end
            if with_envelope
                plot(ax, x, y,'color','black','linewidth',1);
                this.get(fig);
                for i = 1:n
                    %ensure envelope of 95%
                    y_1 = 0*x + true_value + 2*std;
                    y_2 = 0*x + true_value - 2*std;
                    plot(ax, x(1:i), y_1(1:i),'color','black','linewidth',1,'LineStyle','--');
                    plot(ax, x(1:i), y_2(1:i),'color','black','linewidth',1, 'LineStyle','--');
                    this.get(fig);
                end
            end
            this.save(filename=filename,slow=0.01);
            this.exportFrame(filename=filename,ax=ax);
            close(fig);
        end

        function CompareWithExact(this, rots, gt_orientation,filename)
            n = size(gt_orientation,2);
            errors = zeros(n,1);
            t = zeros(n,1);
            dt = 0.05;
            t_tot = 0;
            for i = 1:n
                exact = SO3(quat2rotm(gt_orientation(i)));
                imu_rot = rots{i};
                e = imu_rot - exact;
                errors(i) = norm(e);
                t(i) = t_tot;
                t_tot = t_tot + dt;
            end
            fig = figure;            
            this.plotErrors(t, errors, "t(s)");
            this.exportFrame(filename);
            close(fig);
        end
    end
end