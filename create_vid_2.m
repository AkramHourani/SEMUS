%% create video - Video 2 %%%%%%%%%%%%%%%%%%%%
figure

% 5.a Converting to cartesian coordinates for plotting
[xEast,yNorth,~] = latlon2local(Targetlat,Targetlon,0,GRP);
scatter(xEast(:)/1000,yNorth(:)/1000,2,a(:),'MarkerEdgeColor','none','MarkerFaceColor','flat')

% 5.b Plot geoscatter over the swath

% %% %%create video%%%%%%%%%
% % Create video
VidFilename='[2] Target Points';
v = VideoWriter(VidFilename,'MPEG-4');
v.FrameRate = 15;
v.Quality = 100;
open(v)
% 
% inside a loop
n = 10;
figure
geobasemap satellite
for i = 0:(etaTotal/n)-1 %0:1600/16 

    geoplot(latSwathMidSoI((i*n)+1:((i+1)*n)),lonSwathMidSoI((i*n)+1:((i+1)*n)),'--','LineWidth',1,'MarkerSize',2,'color','w');               % Swath center line (mid swath)
    hold on
    geoplot(GRP(1),GRP(2),'x','LineWidth',1,'MarkerSize',20,'color','r');                          % Swath center point GRP
    geoplot(latSwathL1SoI((i*n)+1:((i+1)*n)),lonSwathL1SoI((i*n)+1:((i+1)*n)),'LineWidth',1.5,'color','w');                                   % Swath edge line 1
    geoplot(latSwathL2SoI((i*n)+1:((i+1)*n)),lonSwathL2SoI((i*n)+1:((i+1)*n)),'LineWidth',1.5,'color','w'); 
    geotickformat -dd

    txt1 = {'- - -','$\times$','------','O'};
    text('Units', 'Normalized', 'Position', [0.77, 0.9, 0], 'string',txt1, 'interpreter','latex','FontSize',10,'BackgroundColor','black','EdgeColor','none','Color','white','horizontalAlignment', 'center' )
    txt2 = {'Swath mid-track','GRP','Swath edges','Testing target points'};
    text('Units', 'Normalized', 'Position', [0.80, 0.9, 0], 'string',txt2, 'interpreter','latex','FontSize',9,'BackgroundColor','black','EdgeColor','none','Color','white','horizontalAlignment', 'left' )

    ax = gca;
    ax.Scalebar.Visible = 'on';
    ax.LatitudeAxis.Label.String='Latitude \circ';
    ax.LongitudeAxis.Label.String='Longitude \circ';
    ax.LatitudeAxis.Label.Color = 'k';
    ax.LongitudeAxis.Label.Color = 'k';
    grid off
    % % title('Satellite swath (target points)','interpreter','latex','FontSize',14);
    gx.Box = 'on';
    set(gca,'LooseInset',get(gca,'TightInset'));

    Frame = getframe(gcf);
    writeVideo(v,Frame);
end
% 
% inside a loop
for i = 1:Param.NtargetsAz
    hold on % Swath edge line 2
    geoscatter(Targetlat(:,i), Targetlon(:,i),'MarkerEdgeColor','y','MarkerFaceColor','none','LineWidth',1)

    Frame = getframe(gcf);
    writeVideo(v,Frame);
end
close (v);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%