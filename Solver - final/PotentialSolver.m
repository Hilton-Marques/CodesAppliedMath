
%
%% Class definition
classdef PotentialSolver < handle
    %%
    % <https://www.mathworks.com/help/matlab/ref/handle-class.html
    % See documentation on *handle* super-class>.
    
    %% Public attributes
    properties (SetAccess = public, GetAccess = public)
        faces = [];
        nfaces = 0;  %Number of faces
        nnode = 0;   % number of nodes
        intnnode = 0; % number of int nodes
        nnel = 0;   % number of elements
        els = [];   % vector of elements
        nodes = []; % vector of nodes
        intnodes = []; % vector of int nodes
        no = 0; % number of points per element
        ndp = 0; % number of prescibred displacement
        ntp = 0; %number of prescribed gradient
        idxDp = [];
        idxTp = [];
        idxUnDp = [];
        idxUnTp = [];
        dp = [];
        tp=[];
        x = [];
        uint = [];
        
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function potential = PotentialSolver(nfaces,no,nnel,nodescor,cord,nelface,nelcord,nodeint,intcord,facescord,q,nodeelement,nnode,elementnode)
            if (nargin > 0)
                potential.nfaces = nfaces;
                potential.no = no;
                faces(1,potential.nfaces) = PotentialFaces();
                for i= 1:potential.nfaces
                    faces(i).k = cond(1,i);
                    faces(i).cord = facescord(i,:);
                    faces(i).q = q(:,i);
                end
                potential.faces = faces;
                potential.nnel = nnel;
                potential.nnode = nnode;
                els(1,potential.nnel) = PotentialElement(); %cria o objeto com parametros nulos
                for i= 1:potential.nnel
                    els(i).face =potential.faces(nelface(1,i));
                    els(i).cord(1,1) = nelcord(1,i);
                    els(i).cord(2,1) = nelcord(2,i);
                    els(i).cord(1,2) = nelcord(3,i);
                    els(i).cord(2,2) = nelcord(4,i);
                    els(i).indnodes = elementnode(:,i); % Dúvida ( como herdar todas as propriedades da classe nós e fornecer a informação de quais nós pertecem ao elemento?
                    els(i).init();
                end
                potential.els = els;
                nodes(1,potential.nnode) = PotentialNode();
                k = 1;
                l=1;
                for i = 1: potential.nnode
                    nodes(i).cord(1) = nodescor(1,i);
                    nodes(i).cord(2) = nodescor(2,i);
                    nodes(i).element = nodeelement(:,i);
                    nodes(i).q(:,2) = potential.els(nodes(i).element(1,1)).q;
                    nodes(i).q(:,1) = potential.els(nodes(i).element(2,1)).q;
                    if nodes(i).q(1,1) == 0 || nodes(i).q(1,2) == 0
                        potential.idxDp(k) = i;
                        if nodes(i).q(1,1) == 0
                            potential.dp(k,1) = nodes(i).q(2,1);
                        else
                            potential.dp(k,1) = nodes(i).q(2,2);
                        end
                        k = k+1;
                    end
                    if i == 1
                        if nodes(i).q(1,2) == 1
                           potential.idxTp(l) = 1;
                            l = l+1;
                        end
                        if nodes(i).q(1,1) == 1
                            potential.idxTp(l)=2*potential.nnode;
                            l = l+1;
                        end
                    else
                         if nodes(i).q(1,1) == 1
                            potential.idxTp(l) = 2*i-1 -1;
                            l = l+1;
                        end
                        if nodes(i).q(1,2) == 1
                            potential.idxTp(l)=2*i - 1;
                            l = l+1;
                        end
                    end
                end
                potential.nodes = nodes;
                potential.ndp = numel(potential.idxDp);
                potential.ntp = numel(potential.idxTp);  
                Q = [potential.nodes.q];
                idxTp = find(Q(1,:)==1);
                potential.tp(:,1) = Q(2,idxTp);
                potential.intnnode = nodeint;
                intnodes(1,potential.nnode) = PotentialNode();
                for i = 1:potential.intnnode
                    intnodes(i).cord(1) = intcord(1,i);
                    intnodes(i).cord(2) = intcord(2,i);
                end
                potential.intnodes = intnodes;
                
            end
            
        end
    end
    methods(Static)
        function circle(x,y,r,c)
            circ = 0 : pi/50 : 2*pi;
            xcirc = x + r * cos(circ);
            ycirc = y + r * sin(circ);
            plot(xcirc, ycirc, 'color', c);
        end
        %Get the shape functions and their derivatives
        function [shapeF,diffShapeF]= shapeFunctions(no)
            if no == 2
                shapeF{1} = @(ksi) 1+(-1).*ksi;
                diffShapeF{1} = @(ksi) (-1);
                shapeF{2} = @(ksi) ksi;
                diffShapeF{2} = @(ksi) (1);
            elseif no ==3
                shapeF{1} = @(ksi) 2.*(-1+ksi).*(-0.5+ksi);
                diffShapeF{1} = @(ksi) -3+4*ksi;
                shapeF{2} = @(ksi) (-4).*(-1+ksi).*ksi;
                diffShapeF{2} = @(ksi) 4-8*ksi;
                shapeF{3} = @(ksi) 2.*(-0.5+ksi).*ksi;
                diffShapeF{3} = @(ksi) -1+4*ksi;
            elseif no ==4
                shapeF{1} = @(ksi) (-1/2).*((-1)+ksi).*((-2)+3.*ksi).*((-1)+3.*ksi);
                diffShapeF{1} = @(ksi) (1/2).*((-11)+36.*ksi+(-27).*ksi.^2);
                shapeF{2} = @(ksi) (9/2).*((-1)+ksi).*ksi.*((-2)+3.*ksi);
                diffShapeF{2} = @(ksi) 9+(-45).*ksi+(81/2).*ksi.^2;
                shapeF{3} = @(ksi) (-9/2).*((-1)+ksi).*ksi.*((-1)+3.*ksi);
                diffShapeF{3} = @(ksi) (-9/2).*(1+(-8).*ksi+9.*ksi.^2);
                shapeF{4} = @(ksi) (1/2).*ksi.*((-2)+3.*ksi).*((-1)+3.*ksi);
                diffShapeF{4} = @(ksi) 1+(-9).*ksi+(27/2).*ksi.^2;
            end
        end
    end
    %% Protected methods
    methods(Access = protected)
        function [xksi,yksi,dxksi,dyksi] = formGeometry(potential,el)
            [shapeF,diffShapeF] = potential.shapeFunctions(potential.no);
            xksi = @(ksi) (0) ;
            yksi = @(ksi) (0) ;
            dxksi = @(ksi) (0) ;
            dyksi = @(ksi) (0) ;
            for k = 1: potential.no
                elementNode = potential.els(el).indnodes(k);
                xksi = @(ksi) xksi(ksi) + potential.nodes(elementNode).cord(1).*shapeF{k}(ksi);
                yksi = @(ksi) yksi(ksi) + potential.nodes(elementNode).cord(2).*shapeF{k}(ksi);
                dxksi = @(ksi) dxksi(ksi) + potential.nodes(elementNode).cord(1).*diffShapeF{k}(ksi);
                dyksi = @(ksi) dyksi(ksi) + potential.nodes(elementNode).cord(2).*diffShapeF{k}(ksi);
            end
        end
        %Calculo a contribuição de G a esquerda e a direita para cada nó
        function GG = formGG(potential,nodes,nnodes)
            [shapeF,~] = PotentialSolver.shapeFunctions(potential.no);
            nnodes(1,nodes+1) = PotentialNode();
            GG =zeros(nodes,2*potential.nnode);
            for i = 1: nodes
                for j = 1:potential.nnode
                    for l =2:-1:1
                        el = potential.nodes(j).element(l);
                        [xksi,yksi,dxksi,dyksi] = formGeometry(potential,el);
                        dr =@(ksi) [dxksi(ksi) ; dyksi(ksi)];
                        r = @(ksi)[xksi(ksi) - nnodes(i).cord(1); yksi(ksi) - nnodes(i).cord(2)];
                        t = (1/(potential.no-1))*(find(potential.els(el).indnodes==j)-1);  %coordenada local do elemento
                        if potential.nodes(j).element(1) ~= potential.nodes(j).element(2)
                            k = [1 potential.no];
                        else
                            pos = find(potential.els(el).indnodes==j);
                            k = [pos pos];
                        end
                        fcn = @(ksi) (vecnorm(dr(t))./(2*pi)).*log(vecnorm(r(ksi)));
                        GG(i,2*j-(l-1)) = integral(@(ksi) fcn(ksi).*shapeF{k(l)}(ksi),0, 1);
                    end
                end
            end
            
            GG = GG(:,[2:2*potential.nnode 1]);
            %conec2 = [potential.nnode*2 1:potential.nnode*2-1];
            
        end
        function HH = formHH(potential,nodes,nnodes)
            [shapeF,~] = potential.shapeFunctions(potential.no);
            indexat = @(expr,index) expr(index);
            nnodes(1,nodes+1) = PotentialNode();  % Why plus one????
            HH=zeros(nodes,potential.nnode);
            for i = 1: nodes
                for j = 1:potential.nnode
                    for l =1:2
                        el = potential.nodes(j).element(l);
                        [xksi,yksi,dxksi,dyksi] = formGeometry(potential,el);
                        r = @(ksi)[xksi(ksi) - nnodes(i).cord(1); yksi(ksi) - nnodes(i).cord(2)];
                        dr = @(ksi) [dxksi(ksi) ; dyksi(ksi)];
                        n = @(ksi) [indexat(dr(ksi),2);-indexat(dr(ksi),1)] ./1;
                        fcn = @(ksi) (( 1 ./ (2*pi*vecnorm(r(ksi)).^2)) .* sum((n(ksi).*r(ksi))) );
                        if potential.nodes(j).element(1) ~= potential.nodes(j).element(2)
                            k = [1 potential.no];
                            HH(i,j) = HH(i,j) + integral(@(ksi) fcn(ksi) .*shapeF{k(l)}(ksi) ,0,1);
                        else
                            pos = find(potential.els(el).indnodes==j); %Encontrando as posições dos nós intermediários
                            k = [pos pos];
                            HH(i,j) = HH(i,j) + 0.5*integral(@(ksi) fcn(ksi) .*shapeF{k(l)}(ksi) ,0,1); % coloquei o 0.5 pois está em um loop
                        end
                    end
                end
            end
        end
        
        function DD = formDD(potential)
            [~,diffShapeF] = PotentialSolver.shapeFunctions(potential.no);
            DD = zeros(potential.nnode,1);
            for i = 1:potential.nnode
                vectors = zeros(2,2);
                for l =1:2
                    el = potential.nodes(i).element(l);
                    dxksi = @(ksi) (0) ;
                    dyksi = @(ksi) (0) ;
                    for k = 1: potential.no
                        dxksi = @(ksi) dxksi(ksi) + potential.nodes(potential.els(el).indnodes(k)).cord(1).*diffShapeF{k}(ksi);
                        dyksi = @(ksi) dyksi(ksi) + potential.nodes(potential.els(el).indnodes(k)).cord(2).*diffShapeF{k}(ksi);
                    end
                    dr = @(ksi) [dxksi(ksi) ; dyksi(ksi)];
                    t = (1/(potential.no-1))*(find(potential.els(el).indnodes==i)-1);
                    v = dr(t)/vecnorm(dr(t));
                    vectors(:,l)=v;
                end
                DD(i,1) = acos(dot(vectors(:,1),-vectors(:,2)))/(2*pi);
            end
        end
    end
    
    %% Public methods
    methods
        function Solve(potential)
            HH = potential.formHH(potential.nnode,potential.nodes);
            DD = potential.formDD();
            WW = DD - sum(HH,2) % Conferência de deslocamento de corpo rígido
            HH = HH - diag(DD);
            GG = potential.formGG(potential.nnode,potential.nodes);
            HHd = zeros(potential.nnode, potential.ndp); 
            GGd = zeros(potential.nnode, potential.ntp);   
            if potential.ndp == potential.nnode
                HHd = HH;
                HHu = [];
            else
                HHd = HH(:,potential.idxDp);
                tot = [1:potential.nnode];
                potential.idxUnDp = setdiff(tot,potential.idxDp,'stable');
                HHu = HH(:,potential.idxUnDp);
            end
            if potential.ntp == 2*potential.nnode
                GGd = GG;
                GGu = [];
            else
                GGd = GG(:,potential.idxTp);
                tot = [1:2*potential.nnode];
                potential.idxUnTp = setdiff(tot,potential.idxTp,'stable');
                GGu = GG(:,potential.idxUnTp);
            end
            A = [HHu -GGu];
            b = -HHd*potential.dp + GGd*potential.tp ;
            potential.x = A\b;
        end
        function Solveint(potential)
            if potential.ndp == potential.nnode
                nd = potential.dp;     
            else
                unkDp = potential.nnode - potential.ndp; % numero de incognitas
                nd([potential.idxUnDp potential.idxDp],1) = [potential.x(1:unkDp); potential.dp];
            end
            if potential.ntp == potential.nnode*2
                nt = potential.tp;
            else
                unkDp = potential.nnode - potential.ndp;
                unkTp = 2*potential.nnode - potential.ntp;
                nt([potential.idxUnTp potential.idxTp],1) = [potential.x(unkDp+1:unkTp+unkDp); potential.tp];
            end
            HH = potential.formHH( potential.intnnode,potential.intnodes);
            GG = potential.formGG( potential.intnnode,potential.intnodes);
            uint = HH*nd - GG*nt
            x=[potential.intnodes.cord];
            uexato = 100.*(1 + x(1:2:end)')'
            erro = abs(uint-uexato)*100./uexato
            potential.uint = uint;
        end
        function DrawMesh(potential)
            % Dúvida : como acessar todos os elementos de um objeto sem
            % precisar do for??
            fac1 = 1 + potential.nnel/40;
            fac2 = 1 + (potential.no-1)/10;
            figure
            hold on
            %x = x(1:2:end);
            for i =1: potential.nfaces
                x(i) = potential.faces(i).cord(1,1);
                y(i) = potential.faces(i).cord(1,2);
            end
            for i = 1:potential.nnode
                xn = potential.nodes(i).cord(1);
                yn = potential.nodes(i).cord(2);
                PotentialSolver.circle(xn,yn,(1/40)/(fac1*fac2),[0 0 1]);
                
                str = sprintf('%.d', i);
                txt = text(xn,yn,str);
                txt.VerticalAlignment = 'top';
                txt.HorizontalAlignment = 'right';
            end
            for i = 1:potential.nnel
                xn = potential.els(i).cord(1,1);
                yn = potential.els(i).cord(2,1);
                xnf = potential.els(i).cord(1,2);
                ynf = potential.els(i).cord(2,2);
                line([xn-((0.03)/fac1), xn+(0.03)/fac1], [(yn-(0.03)/fac1),(yn +(0.03)/fac1)], 'Color', [0 0 0]);
                str = sprintf('%.d', i);
                txt = text((xn+xnf)/2,(yn+ynf)/2,str);
                txt.VerticalAlignment = 'top';
                txt.HorizontalAlignment = 'right';
            end
            h =fill(x,y,[0 0 1]);
            set(h, 'facealpha',0);
            coord = [potential.intnodes.cord];
            hold off
            figure
            hold on
            X = coord(1:2:end)';
            Y = coord(2:2:end)';
            C = potential.uint;
            cmin = min(C);
            cmax = max(C);
            caxis([cmin cmax]);     
            p = patch('Faces',[1 10 20 11 1],'Vertices',[X Y],'FaceVertexCData',C);
p.FaceColor = 'interp';
           cbr = colorbar;
           ylabel(cbr, 'Temperature(u)')
        end
    end
end
