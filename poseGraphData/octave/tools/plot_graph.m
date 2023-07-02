% plot a 2D SLAM graph
function plot_graph(g, iteration)

clf;
hold on;
set(gcf,'color','w');
axis off;

[p, l] = get_poses_landmarks(g);

if (length(l) > 0)
  landmarkIdxX = l;
  landmarkIdxY = l+1;
  plot(g.x(landmarkIdxX), g.x(landmarkIdxY), 'diamond', 'markersize', 8,'markerfacecolor',[0.2588 0.5216 0.9569]);
end

if (length(p) > 0)
  pIdxX = p;
  pIdxY = p+1;
  plot(g.x(pIdxX), g.x(pIdxY), 'o', 'markersize', 8,'markeredge','black','markerfacecolor',"#BDECB6");
end

% draw line segments???
if 1
  poseEdgesP1 = [];
  poseEdgesP2 = [];
  landmarkEdgesP1 = [];
  landmarkEdgesP2 = [];
  for eid = 1:length(g.edges)
    edge = g.edges(eid);
    if (strcmp(edge.type, 'P') ~= 0)
      p1 = g.x(edge.fromIdx:edge.fromIdx+1);
      p2 = g.x(edge.toIdx:edge.toIdx+1);
      poseEdgesP1 = [poseEdgesP1, g.x(edge.fromIdx:edge.fromIdx+1)];
      poseEdgesP2 = [poseEdgesP2, g.x(edge.toIdx:edge.toIdx+1)];
      line([p1(1),p2(1)],[p1(2),p2(2)],'linestyle',':','color','black','linewidth',1);
    elseif (strcmp(edge.type, 'L') ~= 0)
      landmarkEdgesP1 = [landmarkEdgesP1, g.x(edge.fromIdx:edge.fromIdx+1)];
      landmarkEdgesP2 = [landmarkEdgesP2, g.x(edge.toIdx:edge.toIdx+1)];
    end
  end
  
  linespointx = [poseEdgesP1(1,:); poseEdgesP2(1,:)];
  linespointy = [poseEdgesP1(2,:); poseEdgesP2(2,:)];
  
  

  %plot(linespointx, linespointy,'.','color','black');
end

name = "simulation-pose-pose"
exportgraphics(gca,strcat(name,'.png'),'ContentType','vector');

return;
plot(poseEdgesP1(1,:), poseEdgesP1(2,:),'color','r');

%if (columns(poseEdgesP1) > 0)
%end
%if (columns(landmarkEdges) > 0)
%end

hold off;

%figure(1, "visible", "on");
%drawnow;
%pause(0.1);
if (iteration >= 0)
  filename = sprintf('../plots/lsslam_%03d.png', iteration);
  print(filename, '-dpng');
end


end
