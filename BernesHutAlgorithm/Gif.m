classdef Gif < handle
    properties
        filename
        figure
        angleInit_az = 86;
        angleInit_el = 14;
        count_images = 0;
        im = {};
        angle = 2;  
    end
    methods 
        function this = Gif(filename)
%             fig = figure
%             this.figure = fig;
%             this.figure.Position = [100 100 1024  768];
fig = figure( ...
    'color','black', ...
    'menubar','none',...
    'numbertitle', 'off', ...
    'name', 'n-body');
this.figure = fig;
            hold on
            lighting gouraud;
            a1 = camlight('right');
            a2 = camlight('left');
            a3 = camlight('headlight');
            
            set(gca,'XColor', 'none','YColor','none','ZColor','none');
            set(gca,'color','black');
            axis equal;
%             axes( ...
%                 'clipping','off', 'box', 'on',      'ticklength',[0,0], ...
%                 'xticklabel',[],  'yticklabel',[],  'zticklabel',[], ...
%                 'xlim', [-4,4],   'ylim', [-4,4],   'zlim', [-4,4], ...
%                 'xgrid', 'on',    'ygrid', 'on',    'zgrid', 'on', ...
%                 'color',[0,0,0],  'gridcolor',[1,1,0]);
            view(this.angleInit_az,this.angleInit_el)
            this.filename = filename;
        end
        function update(this)
            frame = getframe(this.figure);
            this.im{end+1} = frame2im(frame);
            for i = 1:this.angle
                fRot_az = this.angleInit_az + i/2;
                fRot_el = this.angleInit_el + i/20;
                view(fRot_az,fRot_el);
                %pause(0.05);
                frame = getframe(this.figure);
                this.im{end+1} = frame2im(frame);
            end
            this.angleInit_az = fRot_az;
            this.angleInit_el = fRot_el;
        end
        function save(this)
            for idx = 1:length(this.im)
                [A,map] = rgb2ind(this.im{idx},256);
                if idx == 1
                    imwrite(A,map,this.filename,'gif','LoopCount',Inf,'DelayTime',0.05);
                elseif (idx > 1) && (idx < length(this.im))
                    imwrite(A,map,this.filename,'gif','WriteMode','append','DelayTime',0.0417);
                else
                    imwrite(A,map,this.filename,'gif','WriteMode','append','DelayTime',5);
                end
            end
        end
        function setAngle(this,angle)
            this.angle = angle;
        end
        function exportFrame(~,filename)
            exportgraphics(gca,strcat(filename,'.jpeg'))
        end
    end
end