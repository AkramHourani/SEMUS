%% Create video - Video 3
VidFilename='[3] Unfocused SAR Signal';
v = VideoWriter(VidFilename,'MPEG-4');
v.FrameRate = 25;
v.Quality = 100;
open(v)
% 
% inside a loop
sqd2 = zeros(etaTotal,length(FastTime));% initiate array for draw figure
n = 10;
figure
for i = 0:(etaTotal/n)-1 %0:1600/16 
    sqd2((i*n)+1:((i+1)*n),:) = sqd((i*n)+1:((i+1)*n),:); 
    pc2 = pcolor(FastTime/1e-6,1:etaTotal,abs(sqd2));% sqd = 1600 x 1618
    pc2.LineStyle='none';
    colormap bone
    xlabel('Fast time [\mus]','interpreter','Tex')
    ylabel('Azimuth index','interpreter','Tex')
    gx.Box = 'on';
    set(gca,'LooseInset',get(gca,'TightInset'));
    drawnow

    Frame = getframe(gcf);
    writeVideo(v,Frame);
end
close (v);
hold off
% %%%%%%%%%