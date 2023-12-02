Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);

% 5.a Converting to cartesian coordinates for plotting
[xEast,yNorth,~] = latlon2local(Targetlat,Targetlon,0,GRP);
scatter(xEast(:)/1000,yNorth(:)/1000,2,a(:),'MarkerEdgeColor','none','MarkerFaceColor','flat')
%scatter(xEast(:)/1000,yNorth(:)/1000,20,a(:),'MarkerEdgeColor','k','MarkerFaceColor','none','LineWidth',1)

% 5.b Plot geoscatter over the swath
geoplot(latSwathMidSoI,lonSwathMidSoI,'--','LineWidth',1,'MarkerSize',2,'color','w');               % Swath center line (mid swath)
geobasemap satellite
hold on
geoplot(GRP(1),GRP(2),'x','LineWidth',1,'MarkerSize',20,'color','r');                          % Swath center point GRP
geoplot(latSwathL1SoI,lonSwathL1SoI,'LineWidth',1.5,'color','w');                                   % Swath edge line 1
geoplot(latSwathL2SoI,lonSwathL2SoI,'LineWidth',1.5,'color','w'); 
geotickformat -dd
for i = 1:Param.NtargetsAz
   hold on % Swath edge line 2
   geoscatter(Targetlat(:,i), Targetlon(:,i),'MarkerEdgeColor','y','MarkerFaceColor','none','LineWidth',1)
end

txt1 = {'- - -','$\times$','------','O'};
text('Units', 'Normalized', 'Position', [0.77, 0.9, 0], 'string',txt1, 'interpreter','latex','FontSize',10,'BackgroundColor','black','EdgeColor','none','Color','white','horizontalAlignment', 'center' )
txt2 = {'Swath mid-track','GRP','Swath edges','Testing target points'};
text('Units', 'Normalized', 'Position', [0.80, 0.9, 0], 'string',txt2, 'interpreter','latex','FontSize',9,'BackgroundColor','black','EdgeColor','none','Color','white','horizontalAlignment', 'left' )

ax = gca;
%ax.LatitudeAxis.TickValues=[];
%ax.LongitudeAxis.TickValues=[];
ax.Scalebar.Visible = 'on';
%ax.LatitudeAxis.TickLabels = '';
%ax.LongitudeAxis.TickLabels = '';
%ax.Legend.TextColor = 'white';
%ax.Legend.Color = 'black';
ax.LatitudeAxis.Label.String='Latitude \circ';
ax.LongitudeAxis.Label.String='Longitude \circ';
ax.LatitudeAxis.Label.Color = 'k';
ax.LongitudeAxis.Label.Color = 'k';
grid off
%end here
% 
%comment if overlay
% colormap bone
% axis equal
% hold on
% plot(0,0,'+'); % Mid point (reference)
% xlabel('x-axis [km]','interpreter','latex')
% ylabel('y-axis [km]','interpreter','latex')
% % xticklabels('')
% % yticklabels('')
% % xticks([])
% % yticks([])
% % title('Satellite swath (target points)','interpreter','latex','FontSize',14);
%end here 
gx.Box = 'on';
set(gca,'LooseInset',get(gca,'TightInset'));
drawnow
FilenameG1='Figure6';
print(h_Fig, '-dpng','-r600',FilenameG1)
movefile([FilenameG1 '.png'],'Figures')
