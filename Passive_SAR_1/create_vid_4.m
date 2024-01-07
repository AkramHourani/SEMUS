%% start - create video

% sqd --> unfocused, 1600 x 1618
% src --> range compressed, 1600 x 1618, imagesc(real(src));
% sSLC --> azimuth compressed, 1600 x 1618, J = abs(sSLC);
% imagesc(RangeGround/1000,CrossRange,J)
VidFilename='[4a] Image Processing';
v = VideoWriter(VidFilename,'MPEG-4');
v.FrameRate = 25;
v.Quality = 100;
open(v)
% 
% inside a loop
sqd2 = zeros(etaTotal,length(FastTime));% initiate array for draw figure
src2 = zeros(etaTotal,length(FastTime));
sSLC2 = zeros(etaTotal,length(FastTime));
n = 10;
figure

% from blank to scanning unfocused sar image
% sqd ---> ax.CLim = 1.0e-08 * 0.0001    0.5724
for i = 0:(etaTotal/n)-1 %0:1600/16 

    sqd2((i*n)+1:((i+1)*n),:) = abs(sqd((i*n)+1:((i+1)*n),:)); 

    pc2 = pcolor(FastTime/1e-6,1:etaTotal,sqd2);% sqd = 1600 x 1618
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

VidFilename='[4b] Image Processing';
v = VideoWriter(VidFilename,'MPEG-4');
v.FrameRate = 25;
v.Quality = 100;
open(v)
%from scanning unfocused sar image to range compressed
% src --->  ax.CLim = 1.0e-06 * -0.4886    0.4815
for i = 0:(etaTotal/n)-1 %0:1600/16 

    % method 1 : dynamically changed from unfocused to range compressed
    %src2((i*n)+1:((i+1)*n),:) = real(src((i*n)+1:((i+1)*n),:)); 
    %src2(((i+1)*n)+1:etaTotal,:) = abs(sqd(((i+1)*n)+1:etaTotal,:));

    % method 2 : dynamically changed from blank to range compressed
    src2((i*n)+1:((i+1)*n),:) = real(src((i*n)+1:((i+1)*n),:)); 

    pc2 = pcolor(FastTime/1e-6,1:etaTotal,src2);% sqd = 1600 x 1618
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

VidFilename='[4c] Image Processing';
v = VideoWriter(VidFilename,'MPEG-4');
v.FrameRate = 25;
v.Quality = 100;
open(v)
%from range compressed to azimuth compressed
% sSLC --> ax.Clim = 1.0e-11 * 0.0000    0.3381
for i = 0:(etaTotal/n)-1 %0:1600/16 

    % method 1 : dynamically changed from range compressed to azimuth compressed
    % sSLC2((i*n)+1:((i+1)*n),:) = abs(sSLC((i*n)+1:((i+1)*n),:)); 
    % sSLC2(((i+1)*n)+1:etaTotal,:) = real(src(((i+1)*n)+1:etaTotal,:));

    % method 2 : dynamically changed from blank to azimuth compressed
    sSLC2((i*n)+1:((i+1)*n),:) = abs(sSLC((i*n)+1:((i+1)*n),:)); 

    pc2 = pcolor(FastTime/1e-6,1:etaTotal,sSLC2);% sqd = 1600 x 1618
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
%%end - create video