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

classdef RotationalMechanics < Scene
    properties        
        m_m % prism mass
        m_J % inertia tensor
        m_R % orientation
        m_scale %
        m_w %angular velocity
        m_dt
        m_h
    end
    methods
        function this = RotationalMechanics(l,w,b,m,w0, dt)
            this = this@Scene('rotational.gif','SE3');
            this.m_scale = diag([l,w,b]);
            this.m_m = m;
            this.m_R = SO3;
            this.m_R.setToIdentity();
            this.m_w = w0;
            this.m_dt = dt;
            this.BuildInertialTensor(l,w,b);
            this.Init();
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

        function Solver(this)
            for t = 0:this.m_dt:1
                wd = -inv(this.m_J)*(cross(this.m_w,this.m_J * this.m_w)); % angular acceleration by (3.14)
                this.m_w = this.m_w + wd * this.m_dt;
                this.m_R = this.m_R + this.m_w * this.m_dt; % update
                T = this.m_R.m_data;
                T(4,4) = 1.0;
                delete(this.m_h);
                [this.m_h,pts] = this.drawRobot3D(T,scale=this.m_scale);
                this.get();
            end
            this.save(repeat=true);
        end
    end
end