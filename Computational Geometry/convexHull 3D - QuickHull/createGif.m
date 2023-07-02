classdef createGif < handle
    properties
        filename
        figure
        angleInit_az = 30;
        angleInit_el = -90;
        count_images = 0;
        im = {};
        flag = false;
    end
    methods 
        function this = createGif(fig,filename)
            this.figure = fig;
            view(this.angleInit_az,this.angleInit_el)
            set(gca,'visible','off');
            set(gcf,'color','white');
            this.filename = filename;
        end
        function update(this,rotation)
            if ~this.flag
            frame = getframe(this.figure);
            this.im{end+1} = frame2im(frame);
            for i = 1:2:rotation
                fRot_az = this.angleInit_az + i;
                fRot_el = this.angleInit_el + i/10;
                view(fRot_az,fRot_el);
                pause(0.05);
                frame = getframe(this.figure);
                this.im{end+1} = frame2im(frame);
            end
            this.angleInit_az = fRot_az;
            this.angleInit_el = fRot_el;
            end
        end
        function print(this)
            for idx = 1:length(this.im)
                [A,map] = rgb2ind(this.im{idx},256);
                if idx == 1
                    imwrite(A,map,this.filename,'gif','LoopCount',Inf,'DelayTime',0.01);
                elseif (idx > 1) && (idx < length(this.im))
                    imwrite(A,map,this.filename,'gif','WriteMode','append','DelayTime',0.01);
                end
            end
        end
    end
end