classdef SE2 < Manifold
    properties
    end
    methods
        function this = SE2(data)
            this = this@Manifold();
            if nargin > 0
                this.m_data = data;
                if (size(data,2)) == 1
                    theta = data(3);
                    t1 = data(1);
                    t2 = data(2);
                    this.m_data = [[cos(theta),-sin(theta),t1];[sin(theta),cos(theta),t2];[0,0,1]];                    
                    
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
            inv = [P(1:2,1:2)', -P(1:2,1:2)' * P(1:2,end); P(end,:)];
        end
        function [l, J_l, J_r] = act(this,p)
            P = this.m_data;
            l = P(1:2,3) + P(1:2,1:2)*p;
            [J_l, J_r] = this.getJacobianAct(p);
        end
        function adj = getAdj(this)
            adj = [this.m_data(1:2,1:2),-[-this.m_data(2,3);this.m_data(1,3)];[0,0,1]];
        end
        function adjInv = getAdjInv(this)
            adjInv = [this.m_data(1:2,1:2)',this.m_data(1:2,1:2)'*[-this.m_data(2,3);this.m_data(1,3)];[0,0,1]];
        end
        function [J_l, J_r] = getJacobianAct(this, p)
            J_l = [this.m_data(1:2,1:2), this.m_data(1:2,1:2)*[-p(2);p(1)]];
            J_r = this.m_data(1:2,1:2);
        end
    end
    methods (Static)
        function Q = expMap(v, P)                     
            if nargin == 1
                P = eye(3);
            end
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
        function v = logMap(P,Q)    
            if (nargin == 1)
                C = P;
            else
                invP = [P(1:2,1:2)', -P(1:2,1:2)' * P(1:2,end); P(end,:)];
                C = invP * Q; 
            end            
            %Rotation part
            %Htmp = -logm(H(1:2,1:2)); w = Htmp(1,2);  %implementation more efficient?
            w = atan2(C(2,1), C(1,1));  
            w_sqr = w*w;            
            if (w_sqr <  Manifold.m_epslon)
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
        function P = getRandomElement()            
            v = (rand(2,1)*5 - 2.5);
            w = (rand(1)*2*pi - pi);
            R = [cos(w) -sin(w); sin(w) cos(w)];
            P = [R, v; 0 0 1];
        end
        function v = getRandomTangent(fac)
            if nargin == 0
                fac = 1.0;
            end
            t = fac*(rand(2,1)*2 - 1.0);
            w = fac*(rand(1)*2*pi - pi);
            v = [t;w];
        end
        function Jinv = getJacobianRightInv(v)
            theta = v(3);            
            theta_sqr = theta*theta;
            v1 = v(1);
            v2 = v(2);
            cot_tan_2 = cos(theta/2)/sin(theta/2);
            if (theta_sqr < Manifold.m_epslon)
                a = 2 - theta_sqr/6;
                b = theta/6 + theta^3/360;
            else
                a = theta*cot_tan_2;
                b = (2/theta) - cot_tan_2;
            end
            Jinv = [(0.5)*a,(-0.5)*theta,(0.5)*(v2 + v1*b);...
                    (0.5)*theta,(0.5)*a,(0.5)*(v2*b - v1);0,0,1];
        end
        function v = getTangent(value)

        end
        function n = getTanDim()
            n = 3;
        end
        function n = getDataDim()
            n = [3,3];
        end
        function type = getType()
            type = "SE2";
        end
        function Jleft = getJacobianLeft(v)
        end
        function JinvLeft = getJacobianLeftInv(v)
        end
        function name = getName()
            name = "SE2";
        end
        function dim = getDim()
            dim = 2;
        end
    end
end