Scale = 2.5;
h_Fig = figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2*(2.5/8)/1.618*Scale],'Position',[200 300 800 800*(1.5/8)/1.618*Scale]);

geoplot((SatllaSoI(:,1)),SatllaSoI(:,2),'LineWidth',1.5);hold on  % Satellite subline
geoplot(latSwathMidSoI,lonSwathMidSoI,'--','LineWidth',1,'MarkerSize',2,'color',ColorOrder(5,:));hold on% Swath center line (mid swath)
geoplot(GRP(1),GRP(2),'x','LineWidth',1,'MarkerSize',5,'color',ColorOrder(7,:));hold on   % Swath center point GRP
geoplot(latSwathL1SoI,lonSwathL1SoI,'LineWidth',1.5,'color',ColorOrder(7,:));hold on   % Swath edge line 1
geoplot(latSwathL2SoI,lonSwathL2SoI,'LineWidth',1.5,'color',ColorOrder(7,:));hold on 

geotickformat -dd
%geolimits([-33.55  -33.45],[151 154]) % focus on swath

legend('satellite subtrack','swath mid track','GRP','Swath edges','FontSize',10,'interpreter','latex')
ax=gca;
%ax.LatitudeAxis.TickValues=[];
%ax.LongitudeAxis.TickValues=[];
ax.LatitudeAxis.Label.String='Latitude \circ';
ax.LongitudeAxis.Label.String='Longitude \circ';
ax.Scalebar.Visible = 'on';
%ax.LatitudeAxis.TickLabels = '';
%ax.LongitudeAxis.TickLabels = '';
grid off
drawnow
gx.Box = 'off';
set(gca,'LooseInset',get(gca,'TightInset'));
%title('Swath and Satellite Plot','interpreter','latex','FontSize',14);

FilenameG1='Figure4';
print(h_Fig, '-dpng','-r600',FilenameG1)
movefile([FilenameG1 '.png'],'Figures')