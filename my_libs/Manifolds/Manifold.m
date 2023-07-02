classdef Manifold < handle
    properties (Constant)
        m_alpha = 4; %parameter for schildLadder reduction
        m_epslon = 1e-14;  
    end
    properties       
        m_data = [];
    end
    properties (Abstract)
        %m_data 
    end
    methods
        function this = Manifold()
            %this.m_perm_2 = this.leviCevita(2);
            %this.m_perm_3 = this.leviCevita(3);
        end
        function m = skewG(this, v)
            n = size(v,1);
            if (n == 2)
                m = this.m_perm_2*v;
            elseif (n == 3)
                m = this.m_perm_3*v;
            end
        end
        %schild's Ladder implementation
        %see paper https://hal.inria.fr/hal-02894783 for a upper bound
        %error based on the holonomy 
        function [w,ws] = schildLadder(this,w,P,Q,n_iter)
            if nargin < 5
                n_iter = 20;
            end
            alpha = this.m_alpha;
            t = linspace(0,1,n_iter);
            v = this.logMap(P,Q);
            Qs = this.geodesic(v,P,t);
            ws = zeros(n_iter, size(w,2));
            ws(1,:) = w;
            w_len = norm(w);
            n = n_iter;
            w = w/w_len;
            for i = 2:n_iter
                w = w/(n^alpha);
                P = squeeze(Qs(i-1,:,:));
                Xv = squeeze(Qs(i,:,:));
                Xw = this.expMap(P,w);
                %plot3(pts(:,1),pts(:,2),pts(:,3)); % side
                M = this.expMap(Xw,(0.5)*this.logMap(Xw,Xv));
                %plot3(pts(:,1),pts(:,2),pts(:,3)); % side
                a = this.logMap(P,M);
                Z = this.expMap(P,2*a);
                %plot3(pts(:,1),pts(:,2),pts(:,3)); % side
                w = this.logMap(Xv,Z);
                w = (n^alpha) * w;
                %w = w/norm(w);                
                ws(i,:) = w_len * w;
            end
            w = ws(end,:);
        end
        function pts = geodesic(this,v,P,t)
            if nargin < 4
                t = linspace(0,1,10);
            end
            n = length(t);
            pts = zeros([n,size(P')]);
            for i = 1:n
                ti = t(i);
                pts(i,:,:) = (this.expMap(ti*v,P));
            end
        end
        %P and v define the geodesic, w is the perutbed direction
        %and h is the time
        function j = jacobiField(this, v, P, t, w)
            epslon = this.m_epslon; 
            new_direction = v + epslon*w;
            Q = this.geodesic(v,P,t);
            Q_new = this.geodesic(new_direction,P,t);           
            j = (Q_new - Q) / epslon;
            j_2 = this.logMap(Q_new,Q) / epslon;            
        end
        function w_trans = jacobiAdjoint(this,v,P,t,Y)
            Y_c = this.expMap(t*v,p);
            w = this.logMap(y_c,Y);
            w_trans = this.schildLadder(w,Y_c,P);
        end
        function [x, Y] = samplingPoints(this,P,v,cov,N)
            dim = this.getTanDim();
            cov = cov * eye(dim);
            mean = zeros(1,dim);
            noise = mvnrnd(mean,cov,N);
            x = linspace(0,1,N);
            Y = zeros([N, this.getDataDim()]);
            for i = 1:N
                xi = x(i);
                vi = xi * v;
                Pi = this.expMap(P,vi);
                Y(i,:,:) = this.addNoise(noise(i,:),Pi);
                %Pi = this.exp(vi)*P;
                %Y(i,:) = this.log(this.exp(noise(i,:))*Pi);
            end
        end
        function this = setToIdentity(this)            
            this.m_data = this.getIdentity();
        end
        function I = getIdentity(this)
            v = zeros(this.getTanDim(),1);
            I = this.expMap(v);
        end
        function res = inverse(this)
            res = eval(class(this));
            res.m_data = this.getInverse();
        end
        function res = plus(this, v)
            res = eval(class(this)); 
            res.m_data = this.m_data*this.expMap(v);
        end
        function res = mpower(this,b)
            res = eval(class(this)); 
            res.m_data = this.m_data * b.m_data;
        end
         function t = log(this)
            t = this.logMap(this.m_data);
        end
    end
    methods (Abstract)
        P = addNoise(this,noise,P);
        Y = rep(this,P);
        p = getInverse();
        l = act();
        adj = getAdj(this);
        adjInv = getAdjInv(this);

        J_act = getJacobianAct(this);
    end
    methods (Static, Abstract)
        %fundamental
        Q = expMap(v,P);
        v = logMap(P,Q);
        %utils
        P = getRandomElement();
        v = getRandomTangent(fac);
        n = getTanDim();
        n = getDataDim();
        type = getType();
        Jinv = getJacobianRightInv(v);
        Jleft = getJacobianLeft(v);
        JinvLeft = getJacobianLeftInv(v);
        name = getName();
        dim = getDim();


        %operators
    end 
    methods (Static)        
        function lcMat = leviCevita(N)
            [mats{1:N}] = ndgrid(1:N);
            pairsIndex = nchoosek(1:N,2);
            lcMat = sign(prod(cat(N+1,mats{pairsIndex(:,2)})-...
                  cat(N+1,mats{pairsIndex(:,1)}),N+1));
        end
        function U = skew(v)
            U = [0 -v(3) v(2) ; v(3) 0 -v(1) ; -v(2) v(1) 0 ];
        end
    end
end