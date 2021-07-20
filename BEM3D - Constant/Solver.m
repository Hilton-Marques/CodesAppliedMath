classdef Solver < handle
    properties
        Nodes;
        Elements;
        nNodes;
        nEl;
        ng = 40;
        intNodes;
    end
    methods(Static)
        function [X,Y,Wx,Wy] = triQuad(N)
            v = [0 0; 0 1; 1 0];
            n=1:N;  nnk=2*n+1; A=[1/3 repmat(1,1,N)./(nnk.*(nnk+2))];
            n=2:N; nnk=nnk(n); B1=2/9; nk=n+1; nnk2=nnk.*nnk;
            B=4*(n.*nk).^2./(nnk2.*nnk2-nnk2); ab=[A' [2; B1; B']]; s=sqrt(ab(2:N,2));
            [V,X]=eig(diag(ab(1:N,1),0)+diag(s,-1)+diag(s,1));
            [X,I]=sort(diag(X)); x=(X+1)/2; wx=ab(1,2)*V(1,I)'.^2/4;
            N=N-1; N1=N+1; N2=N+2;  y=cos((2*(N:-1:0)'+1)*pi/(2*N+2));
            L=zeros(N1,N2);  y0=2;  iter=0;
            while max(abs(y-y0))>eps
                L(:,1)=1;    L(:,2)=y;
                for k=2:N1
                    L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
                end
                Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
                y0=y;    y=y0-L(:,N2)./Lp;  iter=iter+1;
            end
            cd=[ 1, 0, 0; -1, 0, 1; 0, 1,-1]*v;
            t=(1+y)/2;  Wx=abs(det(cd(2:3,:)))*wx;  Wy=1./((1-y.^2).*Lp.^2)*(N2/N1)^2;
            [tt,xx]=meshgrid(t,x); yy=tt.*xx;
            X=cd(1,1)+cd(2,1)*xx+cd(3,1)*yy;    Y=cd(1,2)+cd(2,2)*xx+cd(3,2)*yy;
        end
    end
    methods
        function this = Solver(Nodes,Elements,intNodes)
            this.Nodes = Nodes;
            this.nNodes = length(Nodes);
            this.Elements = Elements;
            this.nEl = length(Elements);
            this.intNodes = intNodes; 
        end
        function Calculate(this)
            indexU = [];
            indexQ = [];
            bQ = [];
            bU = [];
            HH = this.FormHH(this.Nodes);
            sum(HH')
            GG = this.FormGG(this.Nodes);
            bU = zeros(this.nEl,1);
            bQ = zeros(this.nEl,1);
            HHtemp = zeros(this.nEl,this.nEl);
            GGtemp = zeros(this.nEl,this.nEl);
            A = zeros(this.nEl,this.nEl);
            for i = 1:this.nEl
                element = this.Elements(i);
                if (element.q(1) == 1)
                    indexU(end+1) = i;
                    bU(i) = element.q(2);
                else
                    indexQ(end+1) = i;
                    bQ(i) = element.q(2);
                end
            end
            if(~isempty(indexQ))
                HHtemp(indexQ,indexQ) = HH(indexQ,indexQ);
            end
            if(~isempty(indexU))
                GGtemp(indexU,indexU) = GG(indexU,indexU);
            end
            b = HH*bU - GG*bQ;
            fileID = fopen('u.txt','w');
            fprintf(fileID,'%12.9f\n',b);
            fclose(fileID);
%             A = GGtemp - HHtemp ;
%             u = A\b
%             nd = zeros(this.nEl,1);
%             nq = zeros(this.nEl,1);
%             nd(indexU,1) = bU(indexU);
%             nd(indexQ,1) = u(indexQ);
%             nq(indexQ,1) = bQ(indexQ);
%             nq(indexU,1) = u(indexU);
%             HH = this.FormHH(this.intNodes); 
%             GG = this.FormGG(this.intNodes);
%             uint = GG*nq - HH*nd
        end
        function HH = FormHH(this,sources)
            nSource = length(sources);
            sources(1,end + 1) = Node();  % Why plus one????
            HH = zeros(nSource,this.nEl);
            [X,Y,Wx,Wy] = this.triQuad(this.ng);
            for i = 1: nSource
                source = sources(i);
                for j = 1:this.nEl
                    field = this.Nodes(j);
                    r = (source.pos - field.pos);
                    d = dot(r,r);
                    if (d > 0)
                        element = field.Element;
                        p1 = element.points(1,:);
                        p2 = element.points(2,:);
                        p3 = element.points(3,:);
                        area = element.area;
                        n = element.n;
                        for l = 1:this.ng
                            for k = 1:this.ng
                                xksi = (1 - X(k,l) - Y(k,l))*p1 + X(k,l)*p2 + Y(k,l)*p3;
                                r = xksi - source.pos;
                                HH(i,j) = HH(i,j) + (area/0.5)*(-dot(r,n)/(4*pi*norm(r)^3))*...
                                    Wx(k)*Wy(l);
                            end
                        end
                    else
                        HH(i,j) = 0.5;
                    end
                end
            end
        end
        function GG = FormGG(this,sources)
            nSource = length(sources);
            sources(1,end + 1) = Node();  % Why plus one????
            GG = zeros(nSource,this.nEl);
            [X,Y,Wx,Wy] = this.triQuad(this.ng);
            for i = 1: nSource
                source = sources(i);
                for j = 1:this.nEl
                    field = this.Nodes(j);
                    element = field.Element;
                    p1 = element.points(1,:);
                    p2 = element.points(2,:);
                    p3 = element.points(3,:);
                    area = element.area;
                    for l = 1:this.ng
                        for k = 1:this.ng
                            xksi = (1 - X(k,l) - Y(k,l))*p1 + X(k,l)*p2 + Y(k,l)*p3;
                            r = xksi - source.pos;
                            GG(i,j) = GG(i,j) + (area/0.5)*(1/(4*pi*norm(r)))*...
                                Wx(k)*Wy(l);
                        end
                    end
                end
            end
        end
    end
    
end