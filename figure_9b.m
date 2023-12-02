Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[1000 150 800 800/1.618*Scale]);

imagesc(RangeGround/1000,CrossRange,J)

ax=gca;
ax.YAxis.Direction = 'reverse';
ax.XAxis.Direction = 'reverse';
xlabel('Ground distance (approx.) [km]')
ylabel('Azimuth [km]')
title('Compressed image')
colormap bone
axis equal
xlim([-1 1]*SwathwidthSoI/2/1000);%% I changed Swathwidth to SwathwidthSoI

% ax = gca;
% %pc = pcolor(Range/1000,etaVec,Img);
% %pc.LineStyle='none';
% ax.YAxis.Direction = 'reverse';
% ax.XAxis.Direction = 'reverse';
% % xticklabels('')
% % yticklabels('')
% xticks('')
% yticks('')
% xlabel('Fast-time','interpreter','Tex')
% ylabel('Azimuth Index','interpreter','Tex')
% title('Final : Compressed Image','interpreter','latex') % Compressed SAR Data, Single Look Complex
% set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',12);
% %colormap turbo
% colormap bone
% % axis equal
% drawnow
%
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',12);
FilenameP6='Figure9b';
print(h_Fig, '-dpng','-r600',FilenameP6)
movefile([FilenameP6 '.png'],'Figures')