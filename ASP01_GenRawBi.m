clc; clear; close all hidden
% This is the main file for simulating space SAR - Passive

%% 0. Load parameters
AP00_ParametersFerdi
Param.NtargetsAz = 11; % number of targets in each eta bin
Param.NtargetsRange = 11; % number of targets in each eta bin

%% 1. Passive SAR
%1a. Create Geometry setup (SoI / Satellite of Interest)
% This Script/function create the satellite orbit
[SatECISoI,SatllaSoI,DateVector] = FP01_CreateSatGeometry(startTime,stopTime,Param,Elem,'Passive_SAR');
etaTotal=length(DateVector); % Total number of slow time steps

%1b. Finding the swath
[latSwathMidSoI,lonSwathMidSoI,slantrangeMidSoI,SwathwidthSoI,latSwathL1SoI,lonSwathL1SoI,latSwathL2SoI,lonSwathL2SoI]=FP02_FindSwath(SatllaSoI,RadPar,E);

%% 2. Active SAR
%2a. Create Geometry setup (Interferer)
startTime = startTime + Delta_t;
stopTime  = startTime + Param.ScanDuration ;
%Elem.TA  = Elem.TA+0.015;  %  in degrees
%Elem.omega  = Elem.omega+0.015;
[SatECII,SatllaI,~] = FP01_CreateSatGeometry(startTime,stopTime,Param,Elem,'Active_SAR'); % This Script/function create the satellite orbit
%SatECII = circshift(SatECISoI,20);
%SatllaI = circshift(SatllaSoI,20);

%2b. Finding the swath of the SoI / Interferer
[latSwathMidI,lonSwathMidI,slantrangeMidI,Swathwidths_mI,latSwathL1I,lonSwathL1I,latSwathL2I,lonSwathL2I] = FP02_FindSwath(SatllaI,RadPar,E);
speed= mean(sqrt(sum((diff(SatECISoI,[],2)).^2)) /Param.ts);
swathlength = speed * time2num(Delta_t);                              % Ground swath length across Azimuth direction

%% 3. Find GRP
% This will find the GRP in the middle of the swath
% find the range migration of the middle of the swath
% This is the index of mid swath
MidEta = round(length(lonSwathL2SoI)/2);
[~,~,RSoI] = geodetic2aer(latSwathMidSoI(MidEta),lonSwathMidSoI(MidEta),0,SatllaSoI(:,1),SatllaSoI(:,2),SatllaSoI(:,3),E);
%RoSoI = min(RSoI);
GRP = [latSwathMidSoI(MidEta),lonSwathMidSoI(MidEta),0]; % ground reference point

% This is the distance between the Active SAR (Int) to the mid of the Passive SAR (SoI) swath
[~,~,RI] = geodetic2aer(latSwathMidSoI(MidEta),lonSwathMidSoI(MidEta),0,SatllaI(:,1),SatllaI(:,2),SatllaI(:,3),E);
%RoI = min(RI);
Ro = min(RSoI+RI);

%% Check the Doppler frequency by checking the maximum velocity of the swath corners
V_max = F10_VelocityCheck(latSwathL1SoI,lonSwathL1SoI,latSwathL2SoI,lonSwathL2SoI,SatllaSoI,E,Ro,Param);

%% Plotting
% Plot the Swath Center line of Passive SAR and Active SAR
% Plot the Swath Edge for Passive SAR and Active SAR
% close all hidden
Scale = 2.5;
figure(1)
h_Fig = figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2*(2.5/8)/1.618*Scale],'Position',[200 300 800 800*(1.5/8)/1.618*Scale]);

geoplot((SatllaSoI(:,1)),SatllaSoI(:,2),'LineWidth',1.5);  % Satellite subline
hold on
geoplot(latSwathMidSoI,lonSwathMidSoI,'--','LineWidth',1,'MarkerSize',2,'color',ColorOrder(5,:));% Swath center line (mid swath)
hold on
geoplot(GRP(1),GRP(2),'x','LineWidth',1,'MarkerSize',5,'color',ColorOrder(7,:));   % Swath center point GRP
hold on
geoplot(latSwathL1SoI,lonSwathL1SoI,'LineWidth',1.5,'color',ColorOrder(7,:));   % Swath edge line 1
hold on
geoplot(latSwathL2SoI,lonSwathL2SoI,'LineWidth',1.5,'color',ColorOrder(7,:)); 
hold on
geotickformat -dd
%geolimits([-33.55  -33.45],[151 154]) % focus on swath

legend('satellite subtrack','swath mid track','GRP','Swath edges','FontSize',10,'interpreter','latex')
ax=gca;
%ax.LatitudeAxis.TickValues=[];
%ax.LongitudeAxis.TickValues=[];
ax.LatitudeAxis.Label.String='Latitude \circ';
ax.LongitudeAxis.Label.String='Longitude \circ';
ax.Scalebar.Visible = 'on';
%ax.LatitudeAxis.TickLabels = '';
%ax.LongitudeAxis.TickLabels = '';
grid off
drawnow
gx.Box = 'off';
set(gca,'LooseInset',get(gca,'TightInset'));
%title('Swath and Satellite Plot','interpreter','latex','FontSize',14);

FilenameG1='Figure4';
print(h_Fig, '-dpng','-r600',FilenameG1)
movefile([FilenameG1 '.png'],'Figures')

%% %%create video - Video 1%%%%%%%%%
% % Create video
VidFilename='[1] Swath and Satellite';
v = VideoWriter(VidFilename,'MPEG-4');
v.FrameRate = 25;
v.Quality = 100;
open(v)
% 
% inside a loop
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
%% %%%%%%%%%%%%%%%%%%%%%%


%% Plot the Swath Edge for Passive SAR and Active SAR
Scale = 1.2;
figure(2)
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);

geoplot(latSwathMidSoI,lonSwathMidSoI,'--','LineWidth',1,'MarkerSize',2,'color','w');               % Swath center line (mid swath)
% geobasemap satellite
hold on
geoplot(GRP(1),GRP(2),'x','LineWidth',1,'MarkerSize',10,'color','w');                          % Swath center point GRP
geoplot(latSwathL1SoI,lonSwathL1SoI,'LineWidth',1.5,'color','w');                                   % Swath edge line 1
geoplot(latSwathL2SoI,lonSwathL2SoI,'LineWidth',1.5,'color','w');                                   % Swath edge line 2
geotickformat -dd
txt1 = {'- - -','$\times$','------'};
text('Units', 'Normalized', 'Position', [0.79, 0.9, 0], 'string',txt1, 'interpreter','latex','FontSize',11,'BackgroundColor','black','EdgeColor','none','Color','white','horizontalAlignment', 'center' )
txt2 = {'Swath mid-track','GRP','Swath edges'};
text('Units', 'Normalized', 'Position', [0.82, 0.9, 0], 'string',txt2, 'interpreter','latex','FontSize',10,'BackgroundColor','black','EdgeColor','none','Color','white','horizontalAlignment', 'left' )
%legend('Swath mid-track','GRP','Swath edges','FontSize',10,'interpreter','latex')
ax = gca;
%ax.LatitudeAxis.TickValues=[];
%ax.LongitudeAxis.TickValues=[];
ax.Scalebar.Visible = 'on';
%ax.LatitudeAxis.TickLabels = '';
%ax.LongitudeAxis.TickLabels = '';
%ax.Legend.TextColor = 'white';
%ax.Legend.Color = 'black';
%ax.AxisColor = [1 1 1];
ax.LatitudeAxis.Label.String='Latitude \circ';
ax.LongitudeAxis.Label.String='Longitude \circ';
ax.LatitudeAxis.Label.Color = 'k';
ax.LongitudeAxis.Label.Color = 'k';
grid off

drawnow
gx.Box = 'on';
set(gca,'LooseInset',get(gca,'TightInset'));
%title('Swath Plot','interpreter','latex','FontSize',14);
FilenameG1='Figure5';
print(h_Fig, '-dpng','-r600',FilenameG1)
movefile([FilenameG1 '.png'],'Figures')
% 
hold on

%% 4. Generate spatial sampling points (Targets)
% Note these are the targets with respect to the SoI (Passive SAR)
[Targetlat,Targetlon]= FP03_GenerateTargets(latSwathL1SoI,lonSwathL1SoI,latSwathL2SoI,lonSwathL2SoI,Param);

%% 5. Get ground reflectivity
a = FP04_GetGroundReflect(Targetlat,Targetlon,latSwathL1SoI,lonSwathL1SoI,latSwathL2SoI,lonSwathL2SoI);

figure(3)
Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);

% 5.a Converting to cartesian coordinates for plotting
[xEast,yNorth,~] = latlon2local(Targetlat,Targetlon,0,GRP);
scatter(xEast(:)/1000,yNorth(:)/1000,2,a(:),'MarkerEdgeColor','none','MarkerFaceColor','flat')
%scatter(xEast(:)/1000,yNorth(:)/1000,20,a(:),'MarkerEdgeColor','k','MarkerFaceColor','none','LineWidth',1)

% 5.b Plot geoscatter over the swath
geoplot(latSwathMidSoI,lonSwathMidSoI,'--','LineWidth',1,'MarkerSize',2,'color','w');               % Swath center line (mid swath)
geobasemap satellite
hold on
geoplot(GRP(1),GRP(2),'x','LineWidth',1,'MarkerSize',20,'color','r');                          % Swath center point GRP
geoplot(latSwathL1SoI,lonSwathL1SoI,'LineWidth',1.5,'color','w');                                   % Swath edge line 1
geoplot(latSwathL2SoI,lonSwathL2SoI,'LineWidth',1.5,'color','w'); 
geotickformat -dd
for i = 1:Param.NtargetsAz
   hold on % Swath edge line 2
   geoscatter(Targetlat(:,i), Targetlon(:,i),'MarkerEdgeColor','y','MarkerFaceColor','none','LineWidth',1)
end

txt1 = {'- - -','$\times$','------','O'};
text('Units', 'Normalized', 'Position', [0.77, 0.9, 0], 'string',txt1, 'interpreter','latex','FontSize',10,'BackgroundColor','black','EdgeColor','none','Color','white','horizontalAlignment', 'center' )
txt2 = {'Swath mid-track','GRP','Swath edges','Testing target points'};
text('Units', 'Normalized', 'Position', [0.80, 0.9, 0], 'string',txt2, 'interpreter','latex','FontSize',9,'BackgroundColor','black','EdgeColor','none','Color','white','horizontalAlignment', 'left' )

ax = gca;
%ax.LatitudeAxis.TickValues=[];
%ax.LongitudeAxis.TickValues=[];
ax.Scalebar.Visible = 'on';
%ax.LatitudeAxis.TickLabels = '';
%ax.LongitudeAxis.TickLabels = '';
%ax.Legend.TextColor = 'white';
%ax.Legend.Color = 'black';
ax.LatitudeAxis.Label.String='Latitude \circ';
ax.LongitudeAxis.Label.String='Longitude \circ';
ax.LatitudeAxis.Label.Color = 'k';
ax.LongitudeAxis.Label.Color = 'k';
grid off
%end here
% 
%comment if overlay
% colormap bone
% axis equal
% hold on
% plot(0,0,'+'); % Mid point (reference)
% xlabel('x-axis [km]','interpreter','latex')
% ylabel('y-axis [km]','interpreter','latex')
% % xticklabels('')
% % yticklabels('')
% % xticks([])
% % yticks([])
% % title('Satellite swath (target points)','interpreter','latex','FontSize',14);
%end here 
gx.Box = 'on';
set(gca,'LooseInset',get(gca,'TightInset'));
drawnow
FilenameG1='Figure6';
print(h_Fig, '-dpng','-r600',FilenameG1)
movefile([FilenameG1 '.png'],'Figures')


% %% create video - Video 2 %%%%%%%%%%%%%%%%%%%%
% figure
% 
% % 5.a Converting to cartesian coordinates for plotting
% [xEast,yNorth,~] = latlon2local(Targetlat,Targetlon,0,GRP);
% scatter(xEast(:)/1000,yNorth(:)/1000,2,a(:),'MarkerEdgeColor','none','MarkerFaceColor','flat')
% 
% % 5.b Plot geoscatter over the swath

% % %% %%create video%%%%%%%%%
% % % Create video
% VidFilename='[2] Target Points';
% v = VideoWriter(VidFilename,'MPEG-4');
% v.FrameRate = 15;
% v.Quality = 100;
% open(v)
% % 
% % inside a loop
% n = 10;
% figure
% geobasemap satellite
% for i = 0:(etaTotal/n)-1 %0:1600/16 
% 
%     geoplot(latSwathMidSoI((i*n)+1:((i+1)*n)),lonSwathMidSoI((i*n)+1:((i+1)*n)),'--','LineWidth',1,'MarkerSize',2,'color','w');               % Swath center line (mid swath)
%     hold on
%     geoplot(GRP(1),GRP(2),'x','LineWidth',1,'MarkerSize',20,'color','r');                          % Swath center point GRP
%     geoplot(latSwathL1SoI((i*n)+1:((i+1)*n)),lonSwathL1SoI((i*n)+1:((i+1)*n)),'LineWidth',1.5,'color','w');                                   % Swath edge line 1
%     geoplot(latSwathL2SoI((i*n)+1:((i+1)*n)),lonSwathL2SoI((i*n)+1:((i+1)*n)),'LineWidth',1.5,'color','w'); 
%     geotickformat -dd
% 
%     txt1 = {'- - -','$\times$','------','O'};
%     text('Units', 'Normalized', 'Position', [0.77, 0.9, 0], 'string',txt1, 'interpreter','latex','FontSize',10,'BackgroundColor','black','EdgeColor','none','Color','white','horizontalAlignment', 'center' )
%     txt2 = {'Swath mid-track','GRP','Swath edges','Testing target points'};
%     text('Units', 'Normalized', 'Position', [0.80, 0.9, 0], 'string',txt2, 'interpreter','latex','FontSize',9,'BackgroundColor','black','EdgeColor','none','Color','white','horizontalAlignment', 'left' )
% 
%     ax = gca;
%     ax.Scalebar.Visible = 'on';
%     ax.LatitudeAxis.Label.String='Latitude \circ';
%     ax.LongitudeAxis.Label.String='Longitude \circ';
%     ax.LatitudeAxis.Label.Color = 'k';
%     ax.LongitudeAxis.Label.Color = 'k';
%     grid off
%     % % title('Satellite swath (target points)','interpreter','latex','FontSize',14);
%     gx.Box = 'on';
%     set(gca,'LooseInset',get(gca,'TightInset'));
% 
%     Frame = getframe(gcf);
%     writeVideo(v,Frame);
% end
% % 
% % inside a loop
% for i = 1:Param.NtargetsAz
%     hold on % Swath edge line 2
%     geoscatter(Targetlat(:,i), Targetlon(:,i),'MarkerEdgeColor','y','MarkerFaceColor','none','LineWidth',1)
% 
%     Frame = getframe(gcf);
%     writeVideo(v,Frame);
% end
% close (v);
% hold off
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% 6. Define The source satellite waveform
%
% [az,elev,slantRange] = geodetic2aer(lat,lon,h,lat0,lon0,h0,spheroid) 
% transforms the geodetic coordinates specified by lat, lon, and h 
% to the local azimuth-elevation-range (AER) spherical coordinates 
% specified by az, elev, and slantRange. 
% Specify the origin of the local AER system with the geodetic coordinates 
% lat0, lon0, and h0. Each coordinate input argument must match the others 
% in size or be scalar. Specify spheroid as the reference spheroid 
% for the geodetic coordinates.
%
%   midEta = 800 
%
%   [~,~,Edge1] = geodetic2aer( -33.8127, 151.3117,0, -33.8906,154.5327,5.1348e+05,E);
%   [~,~,Edge2] = geodetic2aer( -33.8102, 151.2443,0, -33.8906,154.5327,5.1348e+05,E);
%
%   Edge1 = 5.994e+05 ; Edge2 = 6.0279e+05 
%   Swathwidth_SARDistance = 3.3682e+03 m
%   SwathWidthTime = 2.2470e-05 s
%   FastTime = ( -1.3482e-05:1.6667e-08:1.3482e-05) --> 1 x 1618 
%   TimeLength = 1618
%   PulseWidthSamples = 300
%
%
%
%
%
[~,~,Edge1] = geodetic2aer(latSwathL1SoI(MidEta),lonSwathL1SoI(MidEta),0,SatllaSoI(MidEta,1),SatllaSoI(MidEta,2),SatllaSoI(MidEta,3),E);
[~,~,Edge2] = geodetic2aer(latSwathL2SoI(MidEta),lonSwathL2SoI(MidEta),0,SatllaSoI(MidEta,1),SatllaSoI(MidEta,2),SatllaSoI(MidEta,3),E);
% Then calculate Distance and Time
Swathwidth_SARDistance = abs(Edge1-Edge2);
SwathWidthTime = Swathwidth_SARDistance/c*2;
% Calculate FastTime
FastTime = (-SwathWidthTime/2*Param.Margin:RadPar.ts:SwathWidthTime/2*Param.Margin);
TimeLength = length(FastTime);
PulseWidthSamples = round(RadPar.T/(FastTime(end)-FastTime(1))*TimeLength);

%% 7. Generate base chirp (not necessary step, just for testing)
tau = 0;
sb = exp(-1j*pi *   (2*RadPar.fo * tau - RadPar.K*(FastTime-tau).^2   )    ) ...
    .*(FastTime>(-RadPar.T/2+tau)).*(FastTime<(RadPar.T/2+tau));%.*Window; % 1 x 1618

% Plotting the chirp
figure(4)
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

%% 8. (Optional) you can select the Testing value for testing the script
% 0 = optical processing 
% 1 = GRP
% 2 = few targets testing
% 3 = unity reflection
%
% TargetLat --> 11 x 11
% TargetLon --> 11 x 11

Testing = 3 ; 
%FileName = 'SAR_Image1P.mat';

if Testing == 1 % this is for single targets testing
    Targetlat = GRP(1);
    Targetlon = GRP(2);
    a = 1;
    FileName = 'Test01P.mat';
end

NTesting = 5;
if Testing == 2 % this is for three targets testing
    ToPick = randsample(numel(Targetlat),NTesting) ;
    Targetlat = Targetlat(ToPick);
    Targetlon = Targetlon(ToPick);
    a = ones(NTesting,1);
    FileName = 'Test02P.mat';
end

if Testing == 3 % This will force the reflectivity to unity
    a = 1;
    FileName = 'Points.mat';
end

%% 9. Approx azimuth of the satellite to calculate the antenna pattern
%
%   satazSoI = azimuth(-33.8653,154.5331,-33.9160,154.5322,E) +90;
%   satazSoI = 270.7711
%
%   satazI = azimuth(-33.8653,154.5328,-33.9160,154.5320,E) +90;
%   satazI = 270.7711
%
%
if RadPar.Left == 1
    satazSoI = azimuth(SatllaSoI(1,1),SatllaSoI(1,2),SatllaSoI(end,1),SatllaSoI(end,2),E) +90;
    satazI = azimuth(SatllaI(1,1),SatllaI(1,2),SatllaI(end,1),SatllaI(end,2),E) +90;
else
    satazSoI = azimuth(SatllaSoI(1,1),SatllaSoI(1,2),SatllaSoI(end,1),SatllaSoI(end,2),E) -90;
    satazI = azimuth(SatllaI(1,1),SatllaI(1,2),SatllaI(end,1),SatllaI(end,2),E) -90;
end

%% 10. Reference sqd that will be used for template match filtering 
% tauoBi = (RoSoI+RoI)/c;% bi-static delay of the reference point
%
% Ro = 1.2022e+06 m
% tauoBi =  0.0040 s
% reference sqd is at GRP
%
%
%
%
tauoBi  = (Ro)/c;% bi-static delay of the reference point
sqd_ref = zeros(etaTotal,TimeLength);% initialize the reflection matrix, 1600 x 1618
sqd     = (zeros(etaTotal,TimeLength));
disp ('Generating the reference signal...')

for etaIdx=1:etaTotal % 1 : 1600
    sqd_ref(etaIdx,:) = FP05_CalcReflectionBi(1,GRP(1),GRP(2),SatllaI(etaIdx,:),SatllaSoI(etaIdx,:),RadPar,E,mean(satazI),mean(satazSoI),c,tauoBi,FastTime);
end

%% 11. Generate the reflected signal from the entire swath
% This is the longest part of the simulations
% the script will step through the azimuth (slow time) and generate the
% reflected signal from the entire swath
%
%
%
%
%
%
%
tic
disp (['Starting simulation, total steps ',num2str(etaTotal)])
for etaIdx=1:etaTotal
    sqd(etaIdx,:) = FP05_CalcReflectionBi(a,Targetlat,Targetlon,SatllaI(etaIdx,:),SatllaSoI(etaIdx,:),RadPar,E,mean(satazI),mean(satazSoI),c,tauoBi,FastTime);
    disp(etaIdx)
end
toc


%
% sqd --> 1600 x 1618 
%





%%
% close all
% for ctr=1:3
% [~,~,RTargetSoI] = geodetic2aer(Targetlat(ctr),Targetlon(ctr),0,SatllaSoI(:,1),SatllaSoI(:,2),SatllaSoI(:,3),E);
% [~,~,RTargetI] = geodetic2aer(Targetlat(ctr),Targetlon(ctr),0,SatllaI(:,1),SatllaI(:,2),SatllaI(:,3),E);
% ttt = (RTargetI+RTargetSoI)/c-tauoBi;
% plot(ttt)
% hold on
% end

%% 12. Plot the raw unfocused SAR signal (Optional)
figure(5)
%close all hidden
Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
pc = pcolor(FastTime/1e-6,1:etaTotal,abs(sqd));
pc.LineStyle='none';
colormap bone
xlabel('Fast time [\mus]','interpreter','Tex')
ylabel('Azimuth index','interpreter','Tex')
%xticks([])
%yticks([])
%title('Raw time domain (magnitude)','interpreter','latex','FontSize',14);
gx.Box = 'on';
set(gca,'LooseInset',get(gca,'TightInset'));
drawnow
FilenameG4='Figure8';
print(h_Fig, '-dpng','-r600',FilenameG4)
movefile([FilenameG4 '.png'],'Figures')

% %% Create video - Video 3
% VidFilename='[3] Unfocused SAR Signal';
% v = VideoWriter(VidFilename,'MPEG-4');
% v.FrameRate = 25;
% v.Quality = 100;
% open(v)
% % 
% % inside a loop
% sqd2 = zeros(etaTotal,length(FastTime));% initiate array for draw figure
% n = 10;
% figure
% for i = 0:(etaTotal/n)-1 %0:1600/16 
%     sqd2((i*n)+1:((i+1)*n),:) = sqd((i*n)+1:((i+1)*n),:); 
%     pc2 = pcolor(FastTime/1e-6,1:etaTotal,abs(sqd2));% sqd = 1600 x 1618
%     pc2.LineStyle='none';
%     colormap bone
%     xlabel('Fast time [\mus]','interpreter','Tex')
%     ylabel('Azimuth index','interpreter','Tex')
%     gx.Box = 'on';
%     set(gca,'LooseInset',get(gca,'TightInset'));
%     drawnow
% 
%     Frame = getframe(gcf);
%     writeVideo(v,Frame);
% end
% close (v);
% hold off
% % %%%%%%%%%

%%
save(FileName)
%save('Points.mat')
