classdef QuadTree < handle
    properties
        boundary = [];
        capacity = [];
        points = [];
        id = [];
        divided = false;
        level = 0;
        cell_id = 1;
        k = 0;
        northeast = [];
        northwest = [];
        southwest = [];
        southeast = [];
    end
    methods
        function this = QuadTree(boundary,capacity)
            this.boundary = boundary;
            this.capacity = capacity;
        end
        function bool = insert(this,point,id)
            if(~this.boundary.contains(point))
                bool = false;
                return;
            end
            if (size(this.points,1) < this.capacity)
                this.points = [this.points;point];
                this.id = [this.id;id];
                bool = true;
            else
                if(~this.divided)
                    this.subdivide();
                end
                if(this.southwest.insert(point,id))
                    bool = true;
                elseif(this.southeast.insert(point,id))
                    bool = true;
                elseif(this.northwest.insert(point,id))
                    bool = true;
                elseif(this.northeast.insert(point,id))
                    bool = true;
                end
                this.southwest.level = this.level + 1;
                this.southeast.level = this.level + 1;
                this.northwest.level = this.level + 1;
                this.northeast.level = this.level + 1;
            end
        end
        function subdivide(this)
            x = this.boundary.x;
            y = this.boundary.y;
            w = this.boundary.w;
            h = this.boundary.h;
            ne = Rectangle(x + w/2, y + h/2, w/2,h/2);
            this.northeast = QuadTree(ne, this.capacity);
            nw = Rectangle(x - w/2, y + h/2, w/2,h/2);
            this.northwest = QuadTree(nw, this.capacity);
            sw = Rectangle(x - w/2, y - h/2, w/2,h/2);
            this.southwest = QuadTree(sw, this.capacity);
            se = Rectangle(x + w/2, y - h/2, w/2,h/2);
            this.southeast = QuadTree(se, this.capacity);
            this.divided = true;
            %Passing the points to childreen's leaf
            for i = 1:this.capacity
                this.southwest.insert(this.points(i,1:2),this.id(i,1));
                this.southeast.insert(this.points(i,1:2),this.id(i,1));
                this.northwest.insert(this.points(i,1:2),this.id(i,1));
                this.northeast.insert(this.points(i,1:2),this.id(i,1));
            end
        end
        function l_cell_id = label(this,cell_id)
            l_cell_id = label_childreen(this,cell_id);
            this.southwest.label(l_cell_id);
            
            function l_cell_id = label_childreen(this,cell_id)
            if(this.divided)
                if(this.southwest.divided)
                    this.southwest.cell_id = cell_id + 1;
                    cell_id = this.southwest.cell_id;
                end
                if(this.southeast.divided)
                    this.southeast.cell_id = cell_id + 1;
                    cell_id = this.southeast.cell_id;
                end
                if(this.northwest.divided)
                    this.northwest.cell_id = cell_id + 1;
                    cell_id = this.northwest.cell_id;
                end
                if(this.northeast.divided)
                    this.northeast.cell_id = cell_id + 1;
                    cell_id = this.northeast.cell_id;
                end
                l_cell_id = cell_id;
            end
            end
        end
        function query(this)
            
        end
        
        
        function show(this)
            left_corner_x = this.boundary.x - this.boundary.w;
            left_corner_y = this.boundary.y - this.boundary.h;
            rectangle('Position',[left_corner_x,left_corner_y,2*this.boundary.w,...
                2*this.boundary.h]);
            if (this.divided)
                this.northeast.show();
                this.northwest.show();
                this.southwest.show();
                this.southeast.show();
            end
            
            
        end
    end
    
end