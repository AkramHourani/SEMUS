%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main Generation Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Geomtry setup - STEP1.Create Satellite Scenario
[lat, lon]=meshgrid(-37.9292:0.001:-37.8681,144.448:0.01:144.74111);
for i= 1: size(lat,1)
    for j=1:size(lat,2)
    gs(i,j) = groundStation(sc, lat(i,j), lon(i,j));
    end
end
play(sc)
sc.Viewers.CameraReferenceFrame='Inertial';
ac = access(sat, gs(:));
gt = groundTrack(sat);
%% Finding the swath - STEP2.Geometric Simulator
% Plot swath 
% figure(1)
Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);

% geoplot(Satlla(:,1),Satlla(:,2),'LineWidth',1);                                                       % Satellite subline
geoplot((Satlla(:,1)-2.3),Satlla(:,2),'LineWidth',1.5);                                                   % Satellite subline
hold on
geoplot(latSawthMid,lonSwathMid,'--','LineWidth',1,'MarkerSize',2,'color',ColorOrder(5,:));               % Swath center line (mid swath)
geoplot(GRP(1),GRP(2),'x','LineWidth',1,'MarkerSize',5,'color',ColorOrder(7,:));                          % Swath center point GRP
geoplot(latSawthL1,lonSwathL1,'LineWidth',1.5,'color',ColorOrder(7,:));                                   % Swath edge line 1
geoplot(latSawthL2,lonSwathL2,'LineWidth',1.5,'color',ColorOrder(7,:));                                   % Swath edge line 2

% Adding LoRa Transmitter to Geoplot
geoplot(latLORA,lonLORA,'x','LineWidth',1.5,'color',ColorOrder(3,:));                                   % LoRA Transmitter
legend('satellite subtrack','swath mid track','','','','Interferer Tx','FontSize',10,'interpreter','latex')

% legend('satellite subtrack','swath mid track','FontSize',10,'interpreter','latex')
geolimits([-38 -37],[144 146])                                                                            % Latitude - Longitude limits
ax=gca;
ax.LatitudeAxis.TickValues=[];
ax.LongitudeAxis.TickValues=[];
ax.LatitudeAxis.Label.String='Latitude \circ';
ax.LongitudeAxis.Label.String='Longitude \circ';
ax.Scalebar.Visible = 'off';
drawnow
gx.Box = 'on';
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',10);
% title('Swath location','FontSize',14,'interpreter','latex')

% Filename1='Figure8';
% print(h_Fig, '-dpng','-r600',Filename1)
%% Get ground reflectrivity - STEP3.Reflectivity Simulator
% Plot satellite swath
% figure(2)
Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);

pc =pcolor(xEast/1000,yNorth/1000,a);
pc.LineStyle='none';
% scatter(xEast(:)/1000,yNorth(:)/1000,2,a(:),'MarkerEdgeColor','none','MarkerFaceColor','flat')
colormap bone
axis equal
hold on
line(xEast(1,:)/1000,yNorth(1,:)/1000,'LineStyle', '-', 'Color',ColorOrder(7,:), 'LineWidth', 3)
hold on
line(xEast(end,:)/1000,yNorth(end,:)/1000,'LineStyle', '-', 'Color',ColorOrder(7,:), 'LineWidth', 3)
hold on
line(xEast(Idx,:)/1000,yNorth(Idx,:)/1000,'LineStyle', '--', 'Color',ColorOrder(5,:), 'LineWidth', 2)
hold on
plot(0,0,'+','LineWidth',1,'color',ColorOrder(7,:),'MarkerSize', 15);                           % Mid point (reference)
xlabel('x-axis [km]','FontSize',12,'interpreter','latex')
ylabel('y-axis [km]','FontSize',12,'interpreter','latex')
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',12);
% title('Satellite swath (optical)','FontSize',16,'interpreter','latex')

% Filename1='Figure9';
% print(h_Fig, '-dpng','-r600',Filename1)
%% Test antenna pattern - STEP4.Amplitude Simulator
% figure(3)
Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);

% contour(AntennaGain)
imagesc(OffBoreSightAz(:),OffBoreSightRange(:),AntennaGain);
% pc =pcolor(OffBoreSightAz,OffBoreSightRange,AntennaGain);
% pc.LineStyle='none';
% axis equal;
colorbar
xlabel('Azimuth direction [deg]','FontSize',12,'interpreter','latex')
ylabel('Elevation direction [deg]','FontSize',12,'interpreter','latex')
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',12);
% title('Assumed Antenna Gain Pattern','FontSize',16,'interpreter','latex')
% xlim([-1 1])

Filename1='Figure10';
print(h_Fig, '-dpng','-r600',Filename1)
%% Generate base chrip (not nessasry step, just for testing)
figure(4)
plot(FastTime/1e-6,real(sb))
xlabel('Time [\mus]')
ylabel('Real part')
title('Reference pulse [GRP - Mid swath point]')
drawnow
%% STEP5.Waveform Generator
Bar = waitbar( 0 , 'Creating time domain reflections');
for eta=1:etaTotal
    figure(5)
    plot(FastTime/1e-6,real(sqd(eta,:)))
    xlabel('Time [\mus]')
    ylabel('Reflection magnitude')
    disp(eta)
    waitbar(eta/etaTotal,Bar)
end
%% Plotting raw signal
figure(6)
pc =pcolor(FastTime/1e-6,1:etaTotal,abs(sqd));
pc.LineStyle='none';
colormap bone
xlabel('Fast time [ms]')
ylabel('Azimuth index')
title('Raw time domain (magnitude)')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main Processing Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot raw time domain signal
% Define Inverted bone color
% figure(7)
vec = bone;
vec2 = flipud(vec);
Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);

subplot(2,3,1)
pc =pcolor(FastTime/1e-6,1:etaTotal,real(sqd));
pc.LineStyle='none';
colormap(vec2)
% colormap parula
xlabel('Fast time [ms]')
ylabel('Azimuth index','FontSize',8,'interpreter','latex')
title('Step 0:Raw time domain','FontSize',7,'interpreter','latex')
xticks([])
yticks([])
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',10);
%% STEP6.SAR Image Processing
%% Step 1: Range Compression
subplot(2,3,2)
pc =pcolor(FastTime/1e-6,1:etaTotal,real(src));
pc.LineStyle='none';
colormap(vec2)
xlabel('Fast time [ms]')
ylabel('Azimuth index','FontSize',8,'interpreter','latex')
title('Step 1:Range compression','FontSize',7,'interpreter','latex')
xticks([])
yticks([])
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',10);
drawnow
%% Step 2 Azimuth FFT
subplot(2,3,3)
pc =pcolor(FastTime/1e-6,1:etaTotal,abs(S1));
pc.LineStyle='none';
colormap bone
xlabel('Fast time [ms]')
ylabel('Azimuth index','FontSize',8,'interpreter','latex')
title('Step 2:Azimuth FFT','FontSize',7,'interpreter','latex')
xticks([])
yticks([])
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',10);
%% Step 3 Range cell migration compensation
subplot(2,3,4)
plot(1:etaTotal,abs(DeltaR),'Color','k');
xlabel('Azimuth index')
ylabel('Range compensation [m]','FontSize',8,'interpreter','latex')
title('Step 3.1:Range profile','FontSize',7,'interpreter','latex')
xticks([])
yticks([])
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',10);

subplot(2,3,5)
pc =pcolor(FastTime/1e-6,1:etaTotal,real(S2));
pc.LineStyle='none';
colormap(vec2)
xlabel('Fast time [ms]')
ylabel('Azimuth index')
title('Step 3.2:RCMC','FontSize',7,'interpreter','latex')
xticks([])
yticks([])
drawnow
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',10);
%% Step 4 Azimuth compression
subplot(2,3,6)
plot(real(Haz),'Color','k')
xlabel('Azimuth index')
ylabel('Magnitude','FontSize',8,'interpreter','latex')
title('Step 4: Azimuth Matched Filter','FontSize',7,'interpreter','latex')
xticks([])
yticks([])

set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',10);
Filename1='Figure12';
% print(h_Fig, '-dpng','-r600',Filename1)
%% Step 5 Azimuth IFFT
% figure(8) 
RangeBin = 2*RadPar.ts*c;                           % Ground range resolution
Range =(0:RangeBin:(numel(FastTime)-1)*RangeBin)/1000;
speed= mean(sqrt(sum((diff(SatECI,[],2)).^2)) /Param.ts);   % Platform speed = sqrt(Param.mu/(h+Re))
CrossRange = (1:etaTotal)*Param.ts*speed/1000;

Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);

pc =pcolor(Range,CrossRange,abs(sSLC));
% imagesc(abs(sSLC_s{2}))
pc.LineStyle='none';
xlabel('Range [km]')
ylabel('Cross Range [km]')
% title('Step 5: Compressed SAR image using RDA')
% title('Step 5: Focused SAR image with LORA Interference')
colormap bone
% colormap gray
% colormap jet
% xlim([3.9 20.3])
% axis equal
% xlim([0.8 17.05])

set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',10);
% Filename1='Figure17';
% print(h_Fig, '-dpng','-r600',Filename1)
%% Plot interfernce signal
figure(9)
Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
subplot(2,1,1)
spectrogram(signalIQ,500,0,500,fs,'yaxis','centered')
title('LoRa signal injected to the SAR signal');
subplot(2,1,2)
spectrogram(signal_IQ,500,0,500,fs,'yaxis','centered')
title('QPSK modulated signal injected to the SAR signal');
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',10);
Filename1='Figure5';
print(h_Fig, '-dpng','-r600',Filename1)
%% Plot QPSK Signal
plot(t,real(signal_IQ),'r','linewidth',1), grid on;
title('QPSK modulated signal (sum of inphase and Quadrature phase signal)');
xlabel('Time(sec)');
ylabel('Amplitude');
%% Plot Sqd signal with interference
Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
imagesc(real(sqd))
title('Raw Data with LORA Interference')
set(gca,'LooseInset',get(gca,'TightInset'));
%% Plot interference after rearrangement and reshapping
Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
imagesc(real(sLORA))
title('LORA Signal After Rearrangement and Reshapping')
set(gca,'LooseInset',get(gca,'TightInset'));
%% Plot the image the other way around
figure(10)
Range =(0:RangeBin:(numel(FastTime)-1)*RangeBin)/1000;
speed= mean(sqrt(sum((diff(SatECI,[],2)).^2)) /Param.ts);   % Platform speed = sqrt(Param.mu/(h+Re))

% Scale = 1;
% h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);

CrossRange = (1:etaTotal)*Param.ts*speed/1000;
pc =pcolor(CrossRange,Range,abs(sSLC.'));
% imagesc(abs(sSLC))
pc.LineStyle='none';
ylabel('Range [km]')
xlabel('Cross Range [km]')
% title('Step 5: Focused SAR image with LORA Interference')
colormap bone
% colormap gray
% colormap jet
ylim([3.9 20.3])
% axis equal
% xlim([4.3 26])

% set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',10);
% Filename1='Figure17';
% print(h_Fig, '-dpng','-r600',Filename1)
%% Plot stages
subplot(2,3,1);imagesc(real(sqd));title('Raw Data with LORA Interference')
subplot(2,3,2);imagesc(real(src));title('After Range Compression')
subplot(2,3,3);imagesc(real(S1));title('After Azimuth FFT')
subplot(2,3,4);imagesc(real(S2));title('After Range cell migration compensation')
subplot(2,3,5);imagesc(real(S3));title('After Azimuth Compression')
subplot(2,3,6);imagesc(real(sSLC));title('Final Image with LORA Interference')
set(gca,'LooseInset',get(gca,'TightInset'));
sgtitle('Analytical Analysis')
%% Plot Reference
subplot(2,2,1);imagesc(real(sqd_ref));title('Raw Data with LORA Interference')
subplot(2,2,2);imagesc(real(src_ref));title('After Range Compression')
subplot(2,2,3);imagesc(real(S1_ref));title('After Azimuth FFT')
subplot(2,2,4);imagesc(real(S2_ref));title('After Range cell migration compensation')
set(gca,'LooseInset',get(gca,'TightInset'));
sgtitle('Analytical Analysis')
%% Plot stages
subplot(2,3,1);imagesc(abs(sqd));title('Raw Data with LORA Interference')
subplot(2,3,2);imagesc(abs(src));title('After Range Compression')
subplot(2,3,3);imagesc(abs(S1));title('After Azimuth FFT')
subplot(2,3,4);imagesc(abs(S2));title('After Range cell migration compensation')
subplot(2,3,5);imagesc(abs(S3));title('After Azimuth Compression')
subplot(2,3,6);imagesc(abs(sSLC));title('Final Image with LORA Interference')
set(gca,'LooseInset',get(gca,'TightInset'));
sgtitle('Analytical Analysis')
%%
%figure(3)
% subplot(2,1,1)
% clf
% plot(real(S2(:,round(length(S2)/2+2))))
% hold on
% plot(real(Haz)*1e-10)
% plot(imag(sSLC(:,round(length(sSLC)/2+2))))
