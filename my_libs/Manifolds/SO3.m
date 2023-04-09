classdef SO3 < Manifold
    properties
    end
    methods
        function this = SO3()
            this = this@Manifold();
            P = axang2rotm([1,0,0,pi-0.758]);
            w = [0,0,1];
            A = this.schildLadder(w,eye(3),P,40);
            B = P*w';
        end
        function v = addNoise(this,noise_vec,P)
            v = (this.expMap(P,noise_vec));            
        end
        function Y = rep(this, P)
            n = size(P,1);
            Y = zeros(n,3);
            for i = 1:n
                %Y(i,:) = this.logMap(squeeze(P(i,:,:)));
                v = rotm2axang(squeeze(P(i,:,:)));
                Y(i,:) = [v(4)*v(1:3)];
            end
        end
    end
    methods (Static)
        function Q = expMap(v, P)
            if nargin == 1
                P = eye(3);
            end
            U = Manifold.skew(v);
            theta = norm(v);
            theta_sqr = theta*theta;
            if (theta_sqr < Manifold.m_epslon)
                a = 1 - theta_sqr / 6;
                b = 0.5 - theta_sqr / 24;
            else
                a = sin(theta)/theta;
                b = (1 - cos(theta))/theta_sqr;
            end                        
            R = eye(3) + a*U + b*U*U;
            Q = P*R;
        end 
        function Jl = getJacobianLeft(v)
            theta = norm(v);
            theta_sqr = theta*theta;
            theta_cub = theta_sqr*theta;
            if (theta_sqr < Manifold.m_epslon)                
                a = 0.5 - theta_sqr / 24;
                b = 1/6;
            else
                a = (1 - cos(theta))/theta_sqr;
                b = (theta - sin(theta)) / theta_cub;
            end  
            U = Manifold.skew(v);
            Jl = eye(3) + a * U + b * U * U; 
        end       
        function v = logMap(P,Q)
            if (nargin == 1)                
                C = P;
            else
                invP = P';
                C = invP*Q; % P - Q 
            end
            theta = acos((trace(C) - 1)/2);
            if (theta == 0)
                v = [0;0;0];
                return;
            end
            aux = C - C';            
            U = (1/(2*sin(theta)))*aux;
            v = theta*[U(3,2);U(1,3);U(2,1)];
        end
        function JinvL = getJacobianLeftInv(v)
            theta = norm(v);
            theta_sqr = theta*theta;
            if (theta_sqr < Manifold.m_epslon)
                a = (1/12 + theta_sqr/720);
            else
                a = ((1/theta_sqr) - (1 + cos(theta)) / (2*theta*sin(theta)));
            end
            U = Manifold.skew(v);
            JinvL = eye(3) - 0.5 * U + a * U * U;
        end
        function P = getRandomElement()
            P = rand(1,3);
            len = norm(P);
            P = P/len;
            if len > pi
                len = pi - 0.1;
            end
            P = len*P;
            P = axang2rotm([P,len]);
        end
        function v = getRandomVector(P)
            v = rand(1,3);
            len = norm(v);
            v = v/len;
            if len > pi
                len = pi - 0.1;
            end
            v = v*len;
        end
        function n = getTanDim()
            n = 3;
        end
        function n = getDataDim()
            n = [3,3];
        end
        function type = getType()
            type = "SO3";
        end
        
    end
end