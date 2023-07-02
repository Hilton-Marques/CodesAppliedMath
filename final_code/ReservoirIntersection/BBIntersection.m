classdef BBIntersection < handle
    properties
        m_grid
        m_pts
        m_refinement_factor = 2.0;
        m_ncell;
        m_bb;
        m_grid_cell
        m_map
        m_l
        m_collisions
        m_cell_off
        m_ids
    end
    methods
        function this = BBIntersection(pts,ids)
            this.m_ids = ids;
            this.m_grid = pts;
            x = pts(:,1:3:24);
            y = pts(:,2:3:24);
            z = pts(:,3:3:24);
            this.m_pts = [x(:),y(:),z(:)];
            this.m_pts = unique(this.m_pts,'rows');
            this.ComputeRegularDescription();
            this.MapCellToBins();
        end
        function ComputeRegularDescription(this)
            xmin = min(this.m_pts(:,1));
            xmax = max(this.m_pts(:,1));
            ymin = min(this.m_pts(:,2));
            ymax = max(this.m_pts(:,2));
            zmin = min(this.m_pts(:,3));
            zmax = max(this.m_pts(:,3));
            avg_cell_size = this.ComputeAvgCellSize();
            grid_cell = avg_cell_size/this.m_refinement_factor;
            nx = max(1, ceil((xmax - xmin)/grid_cell(1)));
            ny = max(1, ceil((ymax - ymin)/grid_cell(2)));
            nz = max(1, ceil((zmax - zmin)/grid_cell(3)));
            grid_cell = [xmax - xmin, ymax - ymin, zmax - zmin] ./ [nx,ny,nz];
            this.m_grid_cell = grid_cell;
            this.m_ncell = [nx,ny,nz];
            this.m_bb = [[xmin,ymin,zmin];[xmax,ymax,zmax]];
            this.m_l = [xmax - xmin, ymax - ymin,zmax - zmin];
        end
        function avg_side = ComputeAvgCellSize(this)
            x = this.m_grid(:,1:3:24);
            y = this.m_grid(:,2:3:24);
            z = this.m_grid(:,3:3:24);
            x_side = max(x,[],2) - min(x,[],2);
            y_side = max(y,[],2) - min(y,[],2);
            z_side = max(z,[],2) - min(z,[],2);
            x_side_avg = mean(x_side);
            y_side_avg = mean(y_side);
            z_side_avg = mean(z_side);
            avg_side = [x_side_avg,y_side_avg,z_side_avg];            
        end
        function MapCellToBins(this)
            x = this.m_grid(:,1:3:24);
            y = this.m_grid(:,2:3:24);
            z = this.m_grid(:,3:3:24);
            x_side = [min(x,[],2) , max(x,[],2)] - this.m_bb(1,1);
            y_side = [min(y,[],2) , max(y,[],2)] - this.m_bb(1,2);
            z_side = [min(z,[],2) , max(z,[],2)] - this.m_bb(1,3);            
%             i = max(min(floor(x_side / this.m_grid_cell(1)),this.m_ncell(1)-1),0);
%             j = max(min(floor(y_side / this.m_grid_cell(2)),this.m_ncell(2)-1),0);
%             k = max(min(floor(z_side / this.m_grid_cell(3)),this.m_ncell(3)-1),0);
            i = max(min(floor((x_side / this.m_l(1))*this.m_ncell(1)),this.m_ncell(1)-1),0);
            j = max(min(floor((y_side / this.m_l(2))*this.m_ncell(2)),this.m_ncell(2)-1),0);
            k = max(min(floor((z_side / this.m_l(3))*this.m_ncell(3)),this.m_ncell(3)-1),0);
            n = size(i,1);
            map = zeros(n,2,3);
            map(:,:,1) = i;
            map(:,:,2) = j;
            map(:,:,3) = k;   
            this.m_map = map;
            m = prod(this.m_ncell);
            bin_count = zeros(m,1);
            collisions = [];
            cell_ids = [];
            for i = 1:n
                cell = this.GridToBin(this.m_map(i,:,:));
                bin_count(cell) = bin_count(cell) + 1; 
                collisions(end+1:end+length(cell),:) = i; %this.m_ids(i);
                cell_ids(end+1:end+length(cell),:) = cell;
            end
            sum(bin_count);
            bin_offset = zeros(m+1,1);
            count = 1;
            for i = 2:m+1
                count = count + bin_count(i-1);
                bin_offset(i) = count;
            end 
            bins = collisions;
            count = zeros(size(collisions));
            for i = 1:size(collisions,1)
                i;
                cell_id = collisions(i);
                bin_id = cell_ids(i);
                bins(bin_offset(bin_id) + count(bin_id)) = cell_id;
                count(bin_id) = count(bin_id) + 1;
            end
            ids = this.BinToGrid(483004);
            %index = this.GridToBin(ids);
            bb = this.GetBinBB(ids);
            this.m_collisions = bins;
            this.m_cell_off = bin_offset;
        end
        function cells = GridToBin(this, ids)
            ids = reshape(ids,[],3)';
            x = ids(1,1):ids(1,2);
            y = ids(2,1):ids(2,2);
            z = ids(3,1):ids(3,2);
            [X,Y,Z] = meshgrid(x,y,z);
            cells = ((Z * this.m_ncell(2)) + Y)*(this.m_ncell(1)) + X;
            cells = cells(:) + 1;
        end
        function id = BinToGrid(this,index)
            m = this.m_ncell(1);
            n = this.m_ncell(2);
            a = floor(index/(m*n));
            f = index - m*n*a;
            r = mod(f,m);
            i = max(((f - r)/m),1); % + 1
            j = mod(f-1,m); % + 1
            k = a; % + 1
            %id = cat(3,[j,j],[i,i],[k,k]);
            id = [j,i,k];
            check = m*n*(k) + m*(i) + j + 1;
        end
        function [cells,ids] = FindGroup(this,cell)
            bb = this.GetBB(cell);
            ids = this.GetRange(bb)';
            bins = this.GridToBin(ids);  
            x = ids(1,1):ids(2,1);
            y = ids(1,2):ids(2,2);
            z = ids(1,3):ids(2,3);
            [X,Y,Z] = meshgrid(x,y,z);
            ids = [X(:),Y(:),Z(:)];            
            cells = [];
            for i = 1:length(bins)
                if (bins(i) == 1116307)
                    a = 1;
                end
                %check bounding box intersection 
                bb_bin = this.GetBinBB(ids(i,:));
                if ~this.IsBBIntersecting(bb_bin,bb)
                    continue;
                end
                n_cells = this.m_cell_off(bins(i)+1) - this.m_cell_off(bins(i));
                cells(end+1:end+n_cells) = this.m_collisions(this.m_cell_off(bins(i)):this.m_cell_off(bins(i)) + n_cells-1);                
            end
            cells = unique(cells);
        end
        function ids = GetRange(this,bb)
            l_side = bb - this.m_bb(1,:);
            range = (l_side ./ this.m_l) .* this.m_ncell;
            i = max(min(floor(range(:,1)),this.m_ncell(1)-1),0)';
            j = max(min(floor(range(:,2)),this.m_ncell(2)-1),0)';
            k = max(min(floor(range(:,3)),this.m_ncell(3)-1),0)';
            ids = [i;j;k];
        end
        function bb = GetBinBB(this,ids)
            diagonal = this.m_bb(2,:) - this.m_bb(1,:);
            bb_min = this.m_bb(1,:) + ((ids)./this.m_ncell) .* diagonal;
            bb_max = this.m_bb(1,:) + ((ids + 1)./this.m_ncell) .* diagonal;
            bb = [bb_min;bb_max];
        end
    end
    methods (Static)
        function bb = GetBB(pts)
            x = pts(:,1:3:24);
            y = pts(:,2:3:24);
            z = pts(:,3:3:24);
            bb_x = [min(x,[],2) , max(x,[],2)];
            bb_y = [min(y,[],2) , max(y,[],2)];
            bb_z = [min(z,[],2) , max(z,[],2)];
            bb = [bb_x;bb_y;bb_z]';
        end
        function bool = IsBBIntersecting(bb1,bb2)
            bool = false;
            if (bb1(2,1) < bb2(1,1) || bb1(1,1) > bb2(2,1))
                return;
            end
            if (bb1(2,2) < bb2(1,2) || bb1(1,2) > bb2(2,2))
                return;
            end
            if (bb1(2,2) < bb2(1,2) || bb1(1,2) > bb2(2,2))
                return;
            end
            bool = true;
        end
    end
end