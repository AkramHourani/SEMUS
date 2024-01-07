Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);

geoplot(latSwathMidSoI,lonSwathMidSoI,'--','LineWidth',1,'MarkerSize',2,'color','w');               % Swath center line (mid swath)
% geobasemap satellite
hold on
geoplot(GRP(1),GRP(2),'x','LineWidth',1,'MarkerSize',10,'color','w');                          % Swath center point GRP
geoplot(latSwathL1SoI,lonSwathL1SoI,'LineWidth',1.5,'color','w');                                   % Swath edge line 1
geoplot(latSwathL2SoI,lonSwathL2SoI,'LineWidth',1.5,'color','w');                                   % Swath edge line 2
geotickformat -dd
txt1 = {'- - -','$\times$','------'};
text('Units', 'Normalized', 'Position', [0.79, 0.9, 0], 'string',txt1, 'interpreter','latex','FontSize',11,'BackgroundColor','black','EdgeColor','none','Color','white','horizontalAlignment', 'center' )
txt2 = {'Swath mid-track','GRP','Swath edges'};
text('Units', 'Normalized', 'Position', [0.82, 0.9, 0], 'string',txt2, 'interpreter','latex','FontSize',10,'BackgroundColor','black','EdgeColor','none','Color','white','horizontalAlignment', 'left' )
%legend('Swath mid-track','GRP','Swath edges','FontSize',10,'interpreter','latex')
ax = gca;
%ax.LatitudeAxis.TickValues=[];
%ax.LongitudeAxis.TickValues=[];
ax.Scalebar.Visible = 'on';
%ax.LatitudeAxis.TickLabels = '';
%ax.LongitudeAxis.TickLabels = '';
%ax.Legend.TextColor = 'white';
%ax.Legend.Color = 'black';
%ax.AxisColor = [1 1 1];
ax.LatitudeAxis.Label.String='Latitude \circ';
ax.LongitudeAxis.Label.String='Longitude \circ';
ax.LatitudeAxis.Label.Color = 'k';
ax.LongitudeAxis.Label.Color = 'k';
grid off

drawnow
gx.Box = 'on';
set(gca,'LooseInset',get(gca,'TightInset'));
%title('Swath Plot','interpreter','latex','FontSize',14);
FilenameG1='Figure5';
print(h_Fig, '-dpng','-r600',FilenameG1)
movefile([FilenameG1 '.png'],'Figures')
% 