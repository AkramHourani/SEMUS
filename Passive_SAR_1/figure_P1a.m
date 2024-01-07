Scale = 1;
%h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[1000 150 800 800/1.618*Scale]);
h_Fig = figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 2*3.5*2 3.5*2/1.618*Scale],'Position',[1000 150 800 800/1.618*Scale]);
t = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
nexttile
%subplot(2,3,1)
%pc = pcolor(FastTime/1e-6,1:etaTotal,real(sqd));
%pc.LineStyle='none';
imagesc(real(sqd))
colormap bone
xlabel('Fast time [\mus]')
ylabel('Azimuth index')
%xticklabels('')
%yticklabels('')
%xticks([])
%yticks([])

title('(a)', 'Units', 'normalized', 'Position',  [0.05,0.4,0],'interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','none');
%subtitle('Raw time domain (magnitude)')
set(gca,'LooseInset',get(gca,'TightInset'));