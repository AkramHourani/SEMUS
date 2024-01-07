%close all hidden
Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
pc = pcolor(FastTime/1e-6,1:etaTotal,abs(sqd));
pc.LineStyle='none';
colormap bone
xlabel('Fast time [\mus]','interpreter','Tex')
ylabel('Azimuth index','interpreter','Tex')
%xticks([])
%yticks([])
%title('Raw time domain (magnitude)','interpreter','latex','FontSize',14);
gx.Box = 'on';
set(gca,'LooseInset',get(gca,'TightInset'));
drawnow
FilenameG4='Figure8';
print(h_Fig, '-dpng','-r600',FilenameG4)
movefile([FilenameG4 '.png'],'Figures')