%pc = pcolor(FastTime/1e-6,1:etaTotal,abs(S2));
imagesc(abs(S2));
% pc.LineStyle='none';
xlabel('Fast time [\mus]')
ylabel('Azimuth index')
%xticklabels('')
%yticklabels('')
%xticks([])
%yticks([])
title('(d)', 'Units', 'normalized', 'Position',  [0.05,0.4,0],'interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','none');
%subtitle('Step 2: Azimuth FFT')
set(gca,'LooseInset',get(gca,'TightInset'));