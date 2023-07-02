 clear
close(findall(0,'Type','figure'));
clc;

y1 = 1.0e+03 * [3.030313893335954   0.479045796728865   0.474100344843119   0.474099650981131];
y2 = 1.0e+03 * [3.030313893335954   0.477414777912137   0.474104374300813   0.474103823577667];

%show(y1,y2,'simu_pose_land_error');


y1 = 1.0e+06 *[1.795138990772456   0.026716473297156   0.000398317408273   0.000359997114616 0.000359996111601];

y2 = 1.0e+06 * [1.806721036497794   0.004754486376830   0.000360208301020   0.000359972535692 0.000359972530585];

%show(y1,y2,'intel_error');

  y1 = 1.0e+08 * [3.696553355705434   0.639935937145441   1.754502471954035   0.308485117366317   0.004398973482250   0.001137596087356   0.000570428842979   0.000568604377317    0.000568603543455   0.000568603529806]
  
  y2 = 1.0e+10 * [ 0.036965533557054   0.011070586204450   1.154695961657841   0.199072320652164 0.000299690811358   0.000197439328963   0.000009791287576   0.000005720501742 0.000005685960557   0.000005685768074   0.000005685759065   0.000005685758360 0.000005685758298]
  
 %show(y1,y2,'dlr_error');

  
  y1 = 1.0e+08 * [1.388622340753025   0.023216012724990   0.000299862646604   0.000083084946716 0.000082696480638   0.000082694403601   0.000082694243741   0.000082694229105 0.000082694227701];
  
  y2 = 1.0e+08 * [1.462727934716846   0.009834032612710   0.000145552063953   0.000082809011000    0.000082717094248   0.000082717023241   0.000082717022229]
  
   show(y1,y2,'simu_pose_pose_error');
  
  function show(y1,y2,name)
  hold on
  x = linspace(1,size(y1,2),size(y1,2));
  plot(x, y1,'-','color','red','LineWidth',1.0,'marker','*');
  x = linspace(1,size(y2,2),size(y2,2));
  plot(x, y2,'-','color','blue','LineWidth',1.0,'marker','^');
  set(gca,'YScale', 'log') ;
  grid on
  set(gca,'XMinorGrid','off');
  set(gca,'YMinorGrid','off');
  set(gca,'FontSize',12);
  set(gca,'XMinorTick','off');
  set(gca,'YMinorTick','off');
  %set(gca, 'FontName', 'Times New Roman')
  %set(gca,'XTick',x)
  xlabel('$iter$','interpreter','latex');
  ylabel('$\mathbf{F}\left( {\hat{\mathcal{X}}} \right)$','interpreter','latex');
  %ylabel('Error','interpreter','latex');
  axis tight
  legend('linear','manifold','Location','best','interpreter','latex');
  exportgraphics(gca,strcat(name,'.png'),'ContentType','vector');
  end