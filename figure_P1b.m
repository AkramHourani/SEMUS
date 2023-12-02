
plot(real(G));
xlabel('Range Frequency')
ylabel('Real component')
%xticklabels('')
%yticklabels('')
%xticks([])
%yticks([])

title('(b)', 'Units', 'normalized', 'Position',  [0.05,0.4,0],'interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','none');
%subtitle('Real Component of the range matched filter')
set(gca,'LooseInset',get(gca,'TightInset'));
ylim tight
xlim tight
drawnow