plot(1:etaTotal,abs(DeltaR));
xlabel('Azimuth index')
ylabel('Range compensation [m]')
%xticklabels('')
%yticklabels('')
%xticks([])
%yticks([])
title('(e)', 'Units', 'normalized', 'Position',  [0.05,0.4,0],'interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','none');
%subtitle('Step 3.1: Range compensation profile')
set(gca,'LooseInset',get(gca,'TightInset'));