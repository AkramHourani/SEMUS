Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
plot(FastTime/1e-6,real(sb))
xlabel('Time [\mus]','interpreter','Tex')
ylabel('Real part','interpreter','Tex')
xticks([])
yticks([])
title('Reference Pulse [mid swath point]','interpreter','latex','FontSize',12);
gx.Box = 'on';
set(gca,'LooseInset',get(gca,'TightInset'));
drawnow
FilenameG3='FigureG3';
print(h_Fig, '-dpng','-r600',FilenameG3)
movefile([FilenameG3 '.png'],'Figures')
