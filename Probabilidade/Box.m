classdef Box < handle
    properties
        m_p1
        m_p2
        m_name
        m_w
        m_h
        m_color
        m_y_fac
        m_x_fac
        m_area
        m_d
        m_s
        m_m
    end
    methods
        function this = Box(p1,p2,name,color,flag)
            if nargin == 4
                flag = true;
            end
            this.m_p1 = p1;
            this.m_p2 = p2;
            this.m_d = p2 - p1;
            this.m_name = name;
            this.m_w = (p2(1) - p1(1));
            this.m_h = (p2(2) - p1(2));
            this.m_area = this.m_w * this.m_h;
            this.m_color = color;
            this.m_s = this.m_h/this.m_w;
            this.m_m =  this.m_p1 + 0.5*this.m_d;
            if flag
                this.show();
            end
        end
        function [x_text,y_text] = getTextPose(this,caso)
            switch caso
                case 'corner'
                    this.m_y_fac = 0.72;
                    this.m_x_fac = 0.72;
                    %                     if this.m_s > 1
                    %                         m = this.m_h * (1 - this.m_y_fac);
                    %                         this.m_x_fac = 1 - (m/this.m_w);
                    %                     else
                    %                         m = this.m_w * (1 - this.m_x_fac);
                    %                         this.m_y_fac = 1 - (m/this.m_h);
                    %end
                case 'center'
                    this.m_y_fac = 0.5;
                    this.m_x_fac = 0.5;
            end
            x_text = this.m_p1(1) + this.m_x_fac*this.m_d(1);
            y_text = this.m_p1(2) + this.m_y_fac*this.m_d(2);
        end
        function show(this,caso,alpha,fontsize)
            if nargin < 2
                caso = 'corner';
                alpha = '0.4';
                fontsize = 17;
            end
            x = [this.m_p1(1),this.m_p2(1),this.m_p2(1),this.m_p1(1)];
            y = [this.m_p1(2),this.m_p1(2),this.m_p2(2),this.m_p2(2)];
            fill(x,y,this.m_color,'FaceAlpha',alpha);
            [x_text,y_text] = getTextPose(this,caso);
            h = text(x_text,y_text,this.m_name,'FontSize',fontsize,'FontWeight','bold','Interpreter','latex');
            box = h.Extent;
            if strcmp(caso,'corner')
                if (x_text > y_text )
                    m = 1 - x_text;
                    s = box(3) + x_text;
                    off = (s - 1);
                    h.Position = [x_text - off,y_text];
                else
                    m = 1 - y_text;
                    s = box(4) + y_text;
                    off = (s - 1);
                    h.Position = [x_text,y_text-off];
                end
            end
        end
        function Divide(this,name)
            d = this.m_p2 - this.m_p1;
            m = this.m_m;
            p1 = [this.m_p1 + [0,1]*0.5*this.m_d(2)];
            p2 = [this.m_p1 + [1,0]*0.5*this.m_d(1)];
            str = strcat(name,'1}$');
            b1 = Box(p1,p1 + 0.5*d,str,'blue',false);
            str = strcat(name,'2}$');            
            b2 = Box(m,this.m_p2,str,'blue',false);
            str = strcat(name,'3}$');            
            b3 = Box(this.m_p1,m,str,'blue',false);
            str = strcat(name,'4}$');            
            b4 = Box(p2,p2 + 0.5*d,str,'blue',false);
            b1.show('center',0.0,14);
            b2.show('center',0.0,14);
            b3.show('center',0.0,14);
            b4.show('center',0.0,14);
        end
        function res = getCenter(this)
            res = this.m_m;
        end
        function res = getD(this)
            res = norm(this.m_d);
        end
    end
end
