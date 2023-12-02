%pc = pcolor(FastTime/1e-6,1:etaTotal,real(src));
imagesc(real(src))
%pc.LineStyle='none';
xlabel('Fast time [\mus]')
ylabel('Azimuth index')
%xticklabels('')
%yticklabels('')
%xticks([])
%yticks([])

title('(f)', 'Units', 'normalized', 'Position',  [0.05,0.4,0],'interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','none');
%subtitle('Step 3.2: RCMC')
set(gca,'LooseInset',get(gca,'TightInset'));

drawnow
colormap bone
FilenameP5='Figure7';
print(h_Fig, '-dpng','-r600',FilenameP5)
movefile([FilenameP5 '.png'],'Figures')