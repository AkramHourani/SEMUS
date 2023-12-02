Scale = 1.1;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[1000 150 800 800/1.618*Scale]);

clf
% Transform from AzR -> Lat/Lon
etaVec = CalibrationAz:1:(size(Img,1)+CalibrationAz-1);
etaVecM = repmat(etaVec',1,length(Range));
RangeM  = repmat(Range,etaTotal,1);
[LatImg, LonImg] = transformPointsInverse(AzR2LatLon,etaVecM, RangeM);

% Transfrom from Lat/Lon to NEC
[xImg,yImg,~]     = latlon2local(LatImg,LonImg,0,GRP);
%scatter(xImg(:)/1000,yImg(:)/1000,"+","MarkerEdgeColor",ax.ColorOrder(2,:))
hold on

pc = pcolor(xImg/1000,yImg/1000,Img);
%scatter(xEast(:)/1000,yNorth(:)/1000,"+","MarkerEdgeColor",ax.ColorOrder(2,:)) 

pc.LineStyle='none';
colormap turbo
ax.YAxis.Direction = 'reverse';
ax.XAxis.Direction = 'reverse';
%xlim([-1 1]/2* SwathwidthSoI/1000)
axis equal
xlabel('x-axis [km]','interpreter','latex')
ylabel('y-axis [km]','interpreter','latex')
title('Final Step: Projected Image','interpreter','latex')
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',12);
FilenameP7='FigureP7';
print(h_Fig, '-dpng','-r600',FilenameP7)