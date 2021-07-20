
%
%% Class definition
classdef PotentialSolver < handle
    %%
    % <https://www.mathworks.com/help/matlab/ref/handle-class.html
    % See documentation on *handle* super-class>.
    
    %% Public attributes
    properties
        nfaces = 0;  %Number of faces
        nnode = 0;   % number of nodes
        intnnode = 0; % number of int nodes
        nnel = 0;   % number of elements
        els = [];   % vector of elements
        nodes = []; % vector of nodes
        intnodes = []; % vector of int nodes
    end
    properties (Constant)
        ng = 10; % Gauss points
    end
    methods(Static)
        function [x,w]=GaussQuad(N,a,b)
            N = N-1;
            N1 = N+1; N2=N+2;
            xu = linspace(-1,1,N1)';
            % Initial guess
            y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
            % Legendre-Gauss Vandermonde Matrix
            L=zeros(N1,N2);
            % Derivative of LGVM
            Lp=zeros(N1,N2);
            % Compute the zeros of the N+1 Legendre Polynomial
            % using the recursion relation and the Newton-Raphson method
            y0=2;
            % Iterate until new points are uniformly within epsilon of old points
            while max(abs(y-y0))>eps
                L(:,1)=1;
                Lp(:,1)=0;
                L(:,2)=y;
                Lp(:,2)=1;
                for k=2:N1
                    L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
                end
                Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
                y0=y;
                y=y0-L(:,N2)./Lp;
            end
            % Linear map from[-1,1] to [a,b]
            x=(a*(1-y)+b*(1+y))/2;
            % Compute the weights
            w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
        end
        function circle(x,y,r,c)
            circ = 0 : pi/50 : 2*pi;
            xcirc = x + r * cos(circ);
            ycirc = y + r * sin(circ);
            plot(xcirc, ycirc, 'color', c);
        end
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = PotentialSolver(nfaces,no,nnel,nodescor,k,len,nelface,nelcord,nodeint,intcord)
            this.nfaces = nfaces;
            this.nnel = nnel;
            this.nnode = this.nnel*no;
            els(1,this.nnel) = PotentialElement();
            for i= 1:this.nnel
                els(i).k = k(i);
                els(i).len = len(i);
                els(i).nnos = no;
                els(i).face = nelface(1,i);
                els(i).q(1,1) = nelface(2,i);
                els(i).q(2,1) = nelface(3,i);
                els(i).cord(1,1) = nelcord(1,i);
                els(i).cord(2,1) = nelcord(2,i);
                els(i).cord(1,2) = nelcord(3,i);
                els(i).cord(2,2) = nelcord(4,i);
            end
            this.els = els;
            nodes(1,this.nnode) = PotentialNode();
            for i = 1: this.nnode
                nodes(i).cord(1) = nodescor(1,i);
                nodes(i).cord(2) = nodescor(2,i);
            end
            this.nodes = nodes;
            this.intnnode = nodeint;
            intnodes(1,nodeint) = PotentialNode();
            for i = 1:this.intnnode
                intnodes(i).cord(1) = intcord(1,i);
                intnodes(i).cord(2) = intcord(2,i);
            end
            this.intnodes = intnodes;
        end
    end
    
    %% Protected methods
    methods
        function GG = formGG(this,source)
            nSource = length(source);
            source(1,end + 1) = PotentialNode();
            GG = zeros(nSource,this.nnel);
            [ksi,w] = this.GaussQuad(this.ng,0,1);
            for i = 1: nSource
                for j = 1:this.nnode
                    r = (source(i).cord - this.nodes(j).cord);
                    d = dot(r,r);
                    if (d > 0)
                            for k = 1:this.ng
                                xksi = this.els(j).cord(1,1) + (this.els(j).cord(1,2) ...
                                    - this.els(j).cord(1,1))*ksi(k);
                                yksi = this.els(j).cord(2,1) + (this.els(j).cord(2,2) ...
                                    - this.els(j).cord(2,1))*ksi(k);
                                r = [xksi - source(i).cord(1); yksi - source(i).cord(2)];
                                GG(i,j) = GG(i,j) + (this.els(j).len/(2*pi))*log(norm(r))*w(k);
                            end
                    else
                        GG(i,j)= ((this.els(j).len)/(pi*2))*(log(this.els(j).len/2)-1);
                    end
                end
            end
        end
        function HH = formHH(this,source)
            nSource = length(source);
            source(1,end + 1) = PotentialNode();  % Why plus one????
            HH = zeros(nSource,this.nnel);
            [ksi,w] = this.GaussQuad(this.ng,0,1);
            for i = 1: nSource
                for j = 1:this.nnode
                    r = (source(i).cord - this.nodes(j).cord);
                    d = dot(r,r);
                    if (d > 0)
                        for k = 1:this.ng
                            xksi = this.els(j).cord(1,1) + (this.els(j).cord(1,2) ...
                                - this.els(j).cord(1,1))*ksi(k);
                            yksi = this.els(j).cord(2,1) + (this.els(j).cord(2,2) ...
                                - this.els(j).cord(2,1))*ksi(k);
                            r = [xksi - source(i).cord(1); yksi - source(i).cord(2)];
                            n = [(this.els(j).cord(2,2) - this.els(j).cord(2,1)) ; ...
                                -(this.els(j).cord(1,2) - this.els(j).cord(1,1))];
                            n = n/norm(n);
                            HH(i,j) = HH(i,j) + ((this.els(j).len/(2*pi*norm(r)^2))* ...
                                dot(r,n))*w(k);
                        end
                    else
                        HH(i,j)= -0.5;
                    end
                end
            end
            
        end
        function DD = formDD(this)
            DD = zeros(this.nnode,1);
            for i = 1:this.nnel
                vectors = zeros(2,2);
                for l = 1:2
                    if (i == this.nnel && l==2)
                        el = this.els(1);
                        r = el.cord(:,1) - el.cord(:,2);
                        v = r/vecnorm(r);
                        vectors(:,l) = v;
                        continue;
                    end
                    el = this.els(i+l-1);
                    r = el.cord(:,1) - el.cord(:,2);
                    v = r/vecnorm(r);
                    vectors(:,l) = v;
                end
                DD(i,1) = acos(dot(vectors(:,1),-vectors(:,2)))/(2*pi);
            end
        end
    end
    %% Public methods
    methods
        function [conec, dp,tp,x] = Solve(this)
            HH = this.formHH(this.nodes);
            GG = this.formGG(this.nodes);
            k=1;
            g=1;
            HHtot = zeros(this.nnode,this.nnode);
            GGtot = zeros(this.nnode,this.nnode);
            dp = zeros(this.nnode/2,1);
            tp = zeros(this.nnode/2,1);
            conec = zeros(this.nnode/2,2);
            for i = 1: this.nnode
                if this.els(i).q(1,1) == 1 % temperature unknow
                    HHtot(:,k) = HH(:,i);
                    GGtot(:,k+this.nnode/2) = GG(:,i);
                    tp(k,1) = this.els(i).q(2,1);
                    conec(k,1) = i;
                    k = k + 1;
                else %gradient unknown
                    HHtot(:,g+this.nnode/2) = HH(:,i);
                    GGtot(:,g) = GG(:,i);
                    dp(g,1) = this.els(i).q(1,1);
                    conec(g,2) = i;
                    g = g+1;
                end
            end
                        A = [HHtot(:,1:this.nnode/2) -GGtot(:,1:this.nnode/2)];
                        b = -HHtot(:,this.nnode/2+1 :this.nnode)*dp + ...
                            GGtot(:,this.nnode/2+1 :this.nnode)*tp;
%             A = GG;
%             b = HH*dp;
            x = A\b;
        end
        function Solveint(this)
            [conec, dp,tp,x] = Solve(this);
            nd = [ x(1:this.nnode/2,1) ; dp];
            nt = [ x(this.nnode/2 + 1:this.nnode,1); tp];
            nd([conec(:,1) ; conec(:,2)]) = nd;
            nt([conec(:,2) ; conec(:,1)]) = nt;
            HH = this.formHH(this.intnodes);
            GG = this.formGG(this.intnodes);
            uint = HH*nd - GG*nt;
        end
        function DrawMesh(this)
            % Dúvida : como acessar todos os elementos de um objeto sem
            % precisar do for??
            fac1 = 1 + this.nnel/40;
            fac2 = 1 + (1-1)/10;
            
            %x = x(1:2:end);
            %             for i =1: potential.nfaces
            %                 x(i) = potential.faces(i).cord(1,1);
            %                 y(i) = potential.faces(i).cord(1,2);
            %             end
            for i = 1:this.nnode
                xn = this.nodes(i).cord(1);
                yn = this.nodes(i).cord(2);
                PotentialSolver.circle(xn,yn,(1/40)/(fac1*fac2),[0 0 1]);
            end
            for i = 1:this.nnel-1
                xn = this.els(i).cord(1,1);
                yn = this.els(i).cord(2,1);
                xn2 = this.els(i+1).cord(1,1);
                yn2 = this.els(i+1).cord(2,1);
                line([xn-((0.03)/fac1), xn+(0.03)/fac1], [(yn-(0.03)/fac1),(yn +(0.03)/fac1)], 'Color', [0 0 0]);
                line([xn,xn2],[yn,yn2],'Color', [0 0 0]);
                line([this.els(1).cord(1,1),this.els(end).cord(1,1)],[this.els(1).cord(2,1),this.els(end).cord(2,1)],'Color', [0 0 0]);
            end
        end
    end
end
