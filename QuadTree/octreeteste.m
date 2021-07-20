% Clear workspace
clear
close(findall(0,'Type','figure'));

       %pts = rand(10,3)
       pts = [0.4173    0.7803    0.2348; ...
    0.0497    0.3897    0.3532; ...
    0.9027    0.2417    0.8212; ...
    0.9448    0.4039    0.0154; ...
    0.4909    0.0965    0.0430; ...
    0.4893    0.1320    0.1690; ...
    0.3377    0.9421    0.6491; ...
    0.9001    0.9561    0.7317; ...
    0.3692    0.5752    0.6477; ...
    0.1112    0.0598    0.4509];
       OT = OcTree(pts,'binCapacity',5);
       
       OT.shrink
       OT.BinBoundaries
       figure
       boxH = OT.plot;
       cols = lines(OT.BinCount);
       doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
       for i = 1:OT.BinCount
           set(boxH(i),'Color',cols(i,:),'LineWidth', 1+OT.BinDepths(i))
           doplot3(pts(OT.PointBins==i,:),'o','Color',cols(i,:))
       end
       axis image, view(3)

%        points = xlsread('Points','Planilha1');
%        %points = 10*rand(1000,2);
%        draw_obj(points)
%        tic
%        QT = QuadTree2(points,1);
%        QT.cell_count;
%        a = ismember(QT.cell_level,1);
%        i = 1:57;
%        j = i(a);
%        QT.cell_parent
%        toc
%        QT.show()
       
function draw_obj(points)
plot(points(:,1),points(:,2),'o');
for i = 1:size(points,1)
 text(points(i,1),points(i,2),num2str(i),'HorizontalAlignment','left');
end
end
