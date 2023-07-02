classdef SE3 < Manifold
    properties
    end
    methods
        function this = SE3(data)
            this = this@Manifold();
            if nargin > 0
                this.m_data = data;
                [m,n] = size(data);
                if (n == 1 && m == 6)
                    t = data(1:3,1);                    
                    axisang = data(4:end,1);                    
                    this.m_data = [[this.axang2rot(axisang), t];[0,0,0,1]];
                elseif (n == 1 && m == 7)
                    t = data(1:3,1);
                    quat = data(4:end,1);
                    axisang = this.quat2axang(quat);  

                    this.m_data = [[this.axang2rot(axisang), t];[0,0,0,1]];
                    %this.m_data = [[quat2rotm([quat(4);quat(1:3)]'), t];[0,0,0,1]];
                    %this.m_data = [[this.quat2dcm(quat'), t];[0,0,0,1]];
                end
            end

        end
        function v = addNoise(this,noise_vec,P)
            v = (this.expMap(P,noise_vec));            
        end
        function Y = rep(this, P)
            Y = P(:,:,3);
            Y(:,end+1) = acos(P(:,1,1));
        end
        function inv = getInverse(this)
            P = this.m_data;
            inv = [P(1:3,1:3)', -P(1:3,1:3)' * P(1:3,end); P(end,:)];
        end
        function [l, J_l, J_r] = act(this,p)
            P = this.m_data;
            l = P(1:3,4) + P(1:3,1:3)*p;
            [J_l, J_r] = this.getJacobianAct(p);
        end
        function adj = getAdj(this)
            R = this.m_data(1:3, 1:3);
            t = this.m_data(1:3,4);
            adj = [[R,this.skew(t)*R];[0*R, R]];
        end
        function adj = getAdjInv(this)
            R = this.m_data(1:3, 1:3);
            t = this.m_data(1:3,4);
            adj = [[R',-R'*this.skew(t)];[0*R, R']];
        end
        function [J_l, J_r] = getJacobianAct(this, p)
            J_l = [this.m_data(1:3,1:3), -this.m_data(1:3,1:3)*this.skew(p)];
            J_r = this.m_data(1:3,1:3);
        end
    end
    methods (Static)
        function Q = expMap(v, P)                     
            if nargin == 1
                P = eye(4);
            end
            axangle = v(4:end,1);
            t = v(1:3,1);            
            R = SO3.expMap(axangle);
            V = SO3.getJacobianLeft(axangle);
            M = [[R,V*t];[0,0,0,1]];
            Q = P*M;
        end        
        function v = logMap(P,Q)    
            if (nargin == 1)                
                C = P;
            else
                invP = [P(1:3,1:3)', -P(1:3,1:3)' * P(1:3,end); P(end,:)];
                C = invP*Q; % P - Q 
            end
            t = C(1:3,4);
            R = C(1:3,1:3);
            axang = SO3.logMap(R);            
            Jinvl = SO3.getJacobianLeftInv(axang);            
            v = [Jinvl*t;axang];            
        end
        function P = getRandomElement()            
            v = SE3.getRandomTangent(1.0);
            t = v(1:3,1);
            rpy = v(4:6,1);                    
            P = [[SE3.axang2rot(rpy), t];[0,0,0,1]];
        end
        function v = getRandomTangent(fac)
            if nargin == 0
                fac = 1.0;
            end
            t = fac*(rand(3,1)*2 - 1.0);
            w = fac*(rand(1)*2*pi - pi);
            n = rand(3,1);
            n = n/norm(n);
            v = [t;w*n];
        end
        function Jrinv = getJacobianRightInv(v)
            Jrinv = SE3.getJacobianLeftInv(-v);
        end
        function v = getTangent(value)

        end
        function n = getTanDim()
            n = 6;
        end
        function n = getDataDim()
            n = [4,4];
        end
        function type = getType()
            type = "SE2";
        end
        function R = axang2rot(v)
            theta = norm(v);
            if (theta ~= 0)
                v = v/theta;
            end
            U = [0 -v(3) v(2) ; v(3) 0 -v(1) ; -v(2) v(1) 0 ];
            R = eye(3) + sin(theta)*U + (1 - cos(theta))*U*U;
        end
        function v = rot2axang(C)
            theta = acos((trace(C) - 1)/2);
            aux = C - C';
            U = (1/(2*sin(theta)))*aux;
            v = theta*[U(3,2),U(1,3),U(2,1)];
        end
        function Jleft = getJacobianLeft(v)            
            axangle = v(4:6,1);
            Q = SE3.getQ(v);
            Jltheta = SO3.getJacobianLeft(axangle);
            Jleft = [[Jltheta, Q];[0*Jltheta, Jltheta]];
        end
        function Jleftinv = getJacobianLeftInv(v)            
            axangle = v(4:6,1);
            Q = SE3.getQ(v);            
            Jlthetainv = SO3.getJacobianLeftInv(axangle);
            JinvQJinv = Jlthetainv * Q * Jlthetainv; 
            Jleftinv = [[Jlthetainv, -JinvQJinv];[0*Jlthetainv, Jlthetainv]];
        end
        function name = getName()
            name = "SE3";
        end
        function dim = getDim()
            dim = 3;
        end
        function Q = getQ(v)
            t = v(1:3,1);
            axangle = v(4:6,1);
            theta = norm(axangle);
            theta_sqr = theta * theta;
            theta_cub = theta_sqr * theta;
            if theta_sqr < Manifold.m_epslon
                a = (1/6 - theta_sqr/120);
                b = (-1/24 + theta_sqr/720);
                c = (1/120 - theta_sqr/2520);
            else
                a = (theta - sin(theta) )/ theta_cub;
                b = (1 - theta_sqr/2 - cos(theta)) / (theta_sqr * theta_sqr);
                c = -0.5*((1 - theta_sqr/2 - cos(theta))/(theta_sqr*theta_sqr) ...
                    - 3*(theta - sin(theta) - theta_cub/6)/(theta_cub*theta_sqr));
            end
            U = Manifold.skew(t);
            V = Manifold.skew(axangle);
            UV = U*V;
            VU = V*U;
            VUV = V*UV;
            Q = 0.5*U + a*(VU + UV + VUV) - b * (V*VU + UV*V - 3 * VUV)...
                 + c * (VUV*V + V*VUV);
        end
        function v = quat2axang(quat)
            if (quat(4) == 1.0)
                v = [0;0;0];
                return;
            end
            theta = 2 * acos(quat(4));
            den = 1 / sin(theta / 2);
            theta_sqr = theta * theta;
            if (theta_sqr < Manifold.m_epslon)
                den = 2/theta + theta/12;
            end
            x = quat(1) * den;
            y = quat(2) * den;
            z = quat(3) * den;
            v = theta * [x;y;z];
        end
        function dcm = quat2dcm(q)
            qn = bsxfun(@rdivide, q, sqrt(sum(q.^2, 2)));

            dcm = zeros(3,3,size(qn,1));

            dcm(1,1,:) = qn(:,1).^2 + qn(:,2).^2 - qn(:,3).^2 - qn(:,4).^2;
            dcm(1,2,:) = 2.*(qn(:,2).*qn(:,3) + qn(:,1).*qn(:,4));
            dcm(1,3,:) = 2.*(qn(:,2).*qn(:,4) - qn(:,1).*qn(:,3));
            dcm(2,1,:) = 2.*(qn(:,2).*qn(:,3) - qn(:,1).*qn(:,4));
            dcm(2,2,:) = qn(:,1).^2 - qn(:,2).^2 + qn(:,3).^2 - qn(:,4).^2;
            dcm(2,3,:) = 2.*(qn(:,3).*qn(:,4) + qn(:,1).*qn(:,2));
            dcm(3,1,:) = 2.*(qn(:,2).*qn(:,4) + qn(:,1).*qn(:,3));
            dcm(3,2,:) = 2.*(qn(:,3).*qn(:,4) - qn(:,1).*qn(:,2));
            dcm(3,3,:) = qn(:,1).^2 - qn(:,2).^2 - qn(:,3).^2 + qn(:,4).^2;
        end
    end

end