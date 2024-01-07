% %%create video - Video 1%%%%%%%%%
% Create video
VidFilename='[1] Swath and Satellite';
v = VideoWriter(VidFilename,'MPEG-4');
v.FrameRate = 25;
v.Quality = 100;
open(v)

inside a loop
n = 10;
figure
for i = 0:(etaTotal/n)-1 %0:1600/16 
    % geolimits([-34  -33],[151.2 154.6]) % focus on swath
    % hold on
    geoplot((SatllaSoI((i*n)+1:((i+1)*n),1)),SatllaSoI((i*n)+1:((i+1)*n),2),'LineWidth',1.5);  % Satellite subline
    hold on
    geoplot(latSwathMidSoI((i*n)+1:((i+1)*n)),lonSwathMidSoI((i*n)+1:((i+1)*n)),'--','LineWidth',1,'MarkerSize',2,'color',ColorOrder(5,:));% Swath center line (mid swath)
    hold on
    geoplot(GRP(1),GRP(2),'x','LineWidth',1,'MarkerSize',5,'color',ColorOrder(7,:));   % Swath center point GRP
    hold on
    geoplot(latSwathL1SoI((i*n)+1:((i+1)*n)),lonSwathL1SoI((i*n)+1:((i+1)*n)),'LineWidth',1.5,'color',ColorOrder(7,:));   % Swath edge line 1
    hold on
    geoplot(latSwathL2SoI((i*n)+1:((i+1)*n)),lonSwathL2SoI((i*n)+1:((i+1)*n)),'LineWidth',1.5,'color',ColorOrder(7,:)); 
    hold on
    geotickformat -dd
    legend('satellite subtrack','swath mid track','GRP','Swath edges','FontSize',10,'interpreter','latex')
    ax = gca;
    ax.LatitudeAxis.Label.String='Latitude \circ';
    ax.LongitudeAxis.Label.String='Longitude \circ';
    ax.Scalebar.Visible = 'on';
    grid off
    drawnow
    gx.Box = 'off';
    set(gca,'LooseInset',get(gca,'TightInset'));

    Frame = getframe(gcf);
    writeVideo(v,Frame);
end
close (v);
movefile([VidFilename '.mp4'],'Videos')
hold off
% %%%%%%%%%%%%%%%%%%%%%%