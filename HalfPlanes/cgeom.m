classdef cgeom < handle
    properties
    end
    methods
        function this = cgeom()
        end
        function [k2,dualPoints,vol,bool,center] = findInterior(this, M,trans,flag)
            if nargin == 1
                trans = [0,0,0];
                flag = false;
            end
            bool = false;
            %A = A(1:6,:);
            %b = b(1:6,:);
            A = M(:,1:3);
            A(:,4) = 1.0;
            b = -M(:,4);
            n = size(A,2);
            m = size(A,1);
            c = zeros(n,1);
            c(end) = -1;
            tic
            center = linprog(c,A,b);
            r = center(4,1);
            center = center(1:3)';
            plot3(center(1),center(2),center(3),'o','MarkerSize',10,'MarkerFaceColor','green');
            if flag
                x = this.GetChebyshevCenter(M);
            end
            for i = 1:size(M,1)
                n = M(i,1:3);
                d = M(i,4);
                u = center + n * d;
                u = u / norm(u) ;
                angle = dot(u,n);
                if (angle >= 0.0)
                    k2 = 0;
                    dualPoints = [];
                    vol = 0;
                    return;
                end

            end
            bool = true;
            coord_dual = [];
            for i = 1:m
                n = A(i,1:3);
                d = b(i);
                new_d = d - dot(center,n);
                if (new_d == 0)
                    continue;
                end
                coord_dual(end+1,1:3) = ( 1/new_d ) * n;
            end
            [k1,av1] = convhull(coord_dual);
            figure
            hold on
            view(30,30)
            trisurf(k1,coord_dual(:,1),coord_dual(:,2),coord_dual(:,3),'FaceColor','cyan');
            hold off
            volume = this.showDual(coord_dual,k1,center);
            coord_points = coord_dual(unique(k1(:),'stable'),:);
            dualPoints = zeros(size(k1,1),3);
            count = 1;
            for i = 1:size(k1,1)
                inc = k1(i,:);
                p1 = coord_dual(inc(1),:);
                p2 = coord_dual(inc(2),:);
                p3 = coord_dual(inc(3),:);
                normal = cross(p2 - p1, p3 - p1);
                normal = normal / norm(normal);
                d = dot(normal,p1);
                if (d == 0)
                    continue;
                end
                dualPoints(count,:) = (1/d)* normal + center + trans ;
                count = count + 1;
            end
            [k2,vol] = convhull(dualPoints);
            % figure
            % view(30,30);
            % hold on
            %trisurf(k2,dualPoints(:,1),dualPoints(:,2),dualPoints(:,3),'FaceColor','red');
            %hold off;
        end
        function center = GetChebyshevCenter(this,M)
            A = M(:,1:3);
            A(:,4) = 1.0;
            b = -M(:,4);
            n = size(A,2);
            m = size(A,1);
            c = zeros(n,1);
            c(end) = -1;
            newA = A;
            ids = find (b < 0);
            negIds = find(b>=0);
            nIds = length(ids);
            idsM = zeros(size(A,1),nIds);
            for i = 1:nIds
                idsM(ids(i),i) = -1.0;
            end
            newA(ids,:) = -newA(ids,:);
            cAux = zeros(m,1);
            I = eye(m,m);
            newff = [0*c;0*-c;zeros(nIds,1);cAux(negIds);ones(nIds,1)];
            newAA = [newA,-newA,idsM,I(1:m,negIds),-idsM];
            newA = [newA,-newA,idsM,I];
            cAux(ids) = 1.0;
            newf = [0*c;0*-c;zeros(nIds,1);cAux];
            b(ids) = - b(ids);
            n = size(newf,1);
            basis = zeros(1,m);
            count = 0;
            count2 = 1;
            for i = 1 : m
                flag = false;
                for j =  1 : length(ids)
                    if i == ids(j)
                        basis(i) = 20 + j+count;
                        ids(j) = [];
                        flag = true;
                        count = count + 1;
                        break;
                    end
                end
                if flag
                    continue;
                end
                basis(i) = count2 + 8 + nIds;
                count2 = count2 + 1;
            end
            %Part 1
            [x,points,cost,basis] = this.solveLP(newAA,b,newff,basis,b);
            %Part 2
            init = x(basis);
            newf = zeros(20,1);
            n = size(newf,1);
            newAA(:,end-nIds+1:end) = [];
            newf(4,1) = -1.0;
            newf(8,1) = 1.0;
            [x,points,cost,basis] = this.solveLP(newAA,b,newf,basis,init);
            center = x(1:3)' - x(5:7)';
        end

        function [x,points,cost,basis] = solveLP(~,A,b,c,basis,initPoint)
            points = [];
            At = A';
            x = zeros(size(c,1),1);
            d = x;
            dAux = x;
            B = A(:,basis);
            Binv = inv(B);
            x(basis) = initPoint; % xpos at starting corner
            cost = c(basis)'*x(basis);     % cost at starting corner
            n = size(A,2);
            m = size(A,1);
            points(end+1,1:3) = x(1:3,1) - x(4:6,1);
            for iter = 1:100
                d = zeros(size(c,1),1);
                %dAux = zeros(size(c,1),1);
                y = Binv' * c(basis);           % this y may not be feasible
                nonBasics = setdiff(1:n,basis);
                rmin = 0;
                in = -1;
                for k = nonBasics
                    rmin_in = c(k) - At(k,:)*y;
                    if (rmin_in < -.00000001 && rmin_in < rmin)
                        %break
                        rmin = rmin_in;
                        in = k;
                    end
                end
                %[rmin,in] = min(c - A'*y); % minimum r and its index in
                if rmin >= -.00000001      % optimality is reached, r>=0
                    check = A*x - b;
                    break;                 % current x and y are optimal
                end
                Aj = A(:,in);
                d(basis) = Binv * Aj;  % decrease in x from 1 unit of xin
                dAux(basis) =  -d(basis);
                dAux(in) = 1;
                A*dAux;
                xb = x(basis);
                db = d(basis);
                teta = realmax;
                out = -1;
                for i = 1:m
                    dbi = db(i);
                    if (dbi > 1e-3 )
                        inv_dbi = 1 / dbi;
                        xbi = xb(i);
                        value = xbi * inv_dbi;
                        if value < teta
                            teta = value;
                            out = i;
                        end
                    end
                end
                if (out == -1)
                    break;
                end

                cost = cost + teta*rmin  % lower cost at end of step
                x(basis) = x(basis) - teta*d(basis);   % update old x
                x(in) = teta;      % find new positive component of x
                check = A*x - b
                basis(out) = in;      % replace old index by new in basis
                B = A(:,basis);
                BinvL = (1/db(out))*Binv(out,:);
                Q = zeros(12,12);
                BinvBef = Binv;
                for i = 1:12
                    if ( i == out)
                        Binv(i,:) = (1/db(i))*Binv(i,:);
                        %Q(i,i) =  (1/db(i));
                    else
                        Binv(i,:) = Binv(i,:) + (-db(i))*BinvL;
                        Q(i,i) =  1.0;
                        Q(i,out) =  ( -db(i)/db(out) );
                    end
                end
                err = norm(Binv - inv(B));
                points(end+1,1:3) = x(1:3,1) - x(4:6,1);
                Binv = inv(B);
            end
        end
        function volume_other = showDual(~,cell,conec,center)
            volume = 0.0;

            volumes = [];
            volume_other = 0;
            volume_others = [];
            ids = [6,0,5,1,8,2,11];

            ids = ids + 1;

            m = size(cell,1);

            ids_i = 1:m;

            bool = ~ismember(ids_i,ids);

            ordem = [ids,ids_i(bool)];
            %cell = [cell(ids,:);cell(bool,:)];

            order = unique(conec(:),'stable');
            figure
            hold on
            axis equal
            view(30,30)

            %plot3(center(1),center(2),center(3),'o','MarkerSize',10,'MarkerFaceColor','green');
            areas = [];
            for i = order'
                % find star
                star = [];
                for j = 1:size(conec,1)
                    tri = conec(j,:);
                    if ismember(i,tri)
                        star(end+1,1:3) = tri;
                    end
                end
                % put in cylick order
                star = setOrder(star,i);
                % find center of plan
                center = cell(i,:);
                d = norm(center);
                h = 1 / d;
                p0 = (center) * (h * h);
                plot3(p0(1),p0(2),p0(3),'o','MarkerSize',10);
                n = size(star,1)
                %teste
                volume_cell = 0;
                p_first = getDualPoint(cell(star(1,:),:));
                for l = 2:n-1
                    pi = p_first;
                    pj = getDualPoint(cell(star(l,:),:));
                    pk = getDualPoint(cell(star(l+1,:),:));
                    volume_cell = volume_cell + dot(pi, cross(pj,pk))
                    pps = [pi;pj;pk];
                    trisurf([1,2,3],pps(:,1),pps(:,2),pps(:,3));
                end
                volume_others(end+1) = volume_cell;
                volume_other = volume_other + volume_cell;
                area = 0;
                for k = 1:n
                    pi = getDualPoint(cell(star(k,:),:));
                    pj = getDualPoint(cell(star(mod(k,n)+1,:),:));
                    u = pi - p0;
                    v = pj - p0;
                    normal = cross(u,v);
                    area = area + norm(normal)*0.5;
                    pt_triangles = [p0;pi;pj];
                    trisurf([1,2,3],pt_triangles(:,1),pt_triangles(:,2),pt_triangles(:,3));

                end
                areas(end+1) = area;
                volumes(end+1) = area*(h/3);
                volume = volume + area*(h/3);

            end
            hold off
            volume_other = volume_other/6;
        end
        function plane = fittingPlane(~,pts)

            centroide = mean(pts);
            normal = zeros(1,3);
            plane = zeros(1,4);
            %k = convhull(pts(:,1),pts(:,2),pts(:,3));

            %a = getNormal(pts([1,2,3],:));
            % b = getNormal(pts([1,3,4],:));
            % if dot(a,b) < 0
            %     b = -b;
            % end
            % normal = a + b;
            % for i = 1:size(k,2)
            %     normal = normal + getNormal(pts(k(i,:),:));
            % end
            for i = 1:4
                p0 = pts(i,:);
                p1 = pts(mod(i,4) + 1 , :);
                p2 = pts(mod(i+1,4) + 1 , :);
                normal_i = cross(p1-p0,p2-p0);
                len = norm(normal_i);
                if len == 0
                    continue
                end
                normal_i = normal_i / len;
                normal = normal + normal_i;
            end
            if norm(normal) == 0
                return
            end
            normal = normal/norm(normal);
            plane = [normal,-dot(normal, centroide)];
        end
        function out = getDualPoint(~,tri)
            p1 = tri(1,:);
            p2 = tri(2,:);
            p3 = tri(3,:);
            normal = cross(p2 - p1, p3 - p1);
            normal = normal / norm(normal);
            d = dot(normal,p1);
            out = (1/d)* normal ;
        end
    end
end