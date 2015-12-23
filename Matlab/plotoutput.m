% Plot

% subplot(2,2,1);
subplot(1,2,1);
plot(time,stran(3,:),'-b','LineWidth',4);
xlabel('Time [s]');
ylabel('Strain component 33 [-]');
title('Strain vs. Time');
set (gcf,'windowstyle','normal');  %gcf-figure gca-axis
set (gcf,'Position',[500,300,400,260]);   %left bottom weight height
set(gca,'FontName','Times New Roman','FontSize',27,'LineWidth',1.0);
set(gca,'Position',get(gca,'Position')+[0 0.02 0 0])

% subplot(2,2,2);
subplot(1,2,2);
plot(stran(3,:),stress(3,:),'-b', eps33_an, S33_an,'xr', 'LineWidth',4,'MarkerSize',10); 
xlabel('Strain component 33 [-]');
ylabel('Stress component 33 [MPa]');
title('Stress vs. Strain');
set (gcf,'windowstyle','normal');  %gcf-figure gca-axis
set (gcf,'Position',[500,300,400,260]);   %left bottom weight height
set(gca,'FontName','Times New Roman','FontSize',27,'LineWidth',1.0);
set(gca,'Position',get(gca,'Position')+[0 0.02 0 0])

% subplot(2,2,3);
% % subplot(1,2,1);
% plot(time,stran(6,:),'-b','LineWidth',4);
% xlabel('Time [s]');
% ylabel('Strain component 23 [-]');
% title('Strain vs. Time');
% set (gcf,'windowstyle','normal');  %gcf-figure gca-axis
% set (gcf,'Position',[500,300,400,260]);   %left bottom weight height
% set(gca,'FontName','Times New Roman','FontSize',27,'LineWidth',1.0);
% set(gca,'Position',get(gca,'Position')+[0 0.02 0 0])
% 
% subplot(2,2,4);
% % subplot(1,2,2);
% plot(stran(6,:),stress(6,:),'-b', gam23_an, S23_an,'xr','LineWidth',4,'MarkerSize',10);
% xlabel('Strain component 23 [-]');
% ylabel('Stress component 23 [MPa]');
% title('Stress vs. Strain');
% set (gcf,'windowstyle','normal');  %gcf-figure gca-axis
% set (gcf,'Position',[500,300,400,260]);   %left bottom weight height
% set(gca,'FontName','Times New Roman','FontSize',27,'LineWidth',1.0);
% set(gca,'Position',get(gca,'Position')+[0 0.02 0 0])