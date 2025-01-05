clc; clear; close all
close all hidden;
load('SAR_Image')                     % That is my final image with azimuth 5°
%% This is a raw-wise FFT / IFFT
fft1d2 = @ (x) fftshift(fft(fftshift(x,2),[],2),2);
ifft1d2 = @ (x) ifftshift(ifft(ifftshift(x,2),[],2),2);
% This is a cloumn-wise FFT - Azimuth
fft1d1 = @ (x) fftshift(fft(fftshift(x,1),[],1),1);
ifft1d1 = @ (x) ifftshift(ifft(ifftshift(x,1),[],1),1);
%% Add noise and interference to the received signal
A02_Parameters                                       % Load interference parameters
% load('sLORA.mat')
% % Add AWGN to the recieved signal
NI01_GenerateAWGN
sqd = sqd + AWGN;                                    % Signal to interference = Noise.SNR
% Add LORA signal
% NI02_GenerateLORA
% sqd = sqd + sLORA;                                   % Signal to interference = LORA.SIR
% Add AM signal
% NI03_GenerateAM
% sqd = sqd + sAM;                                    % Signal to interference = AM.SIR
% % Add QPSK signal
% NI04_GenerateQPSK
% sqd = sqd + sQPSK;                                  % Signal to interference = QPSK.SIR
% Add Radar signal
% NI05_GenerateRadarTx
% sqd = sqd + sInfR;                                   % Signal to interference = IR.SIR
% figure,imagesc(real(sqd))
sqdT = sqd.';
figure;pwelch(sqdT(:),size(sqd,1),[],size(sqd,1),RadPar.fs,'centered')
%% Plot Optical reflectrivity Figure210
Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
latEdge = [min([latSwathL1;latSwathL2]), max([latSwathL1;latSwathL2])];
lonEdge = [min([lonSwathL1;lonSwathL2]), max([lonSwathL1;lonSwathL2])];
geolimits(latEdge,lonEdge)
% geolimits([-36.1664  -36.0700],[150.0794  150.1037])

geobasemap satellite
ax=gca;
% pc.LineStyle='none';
grid on
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',12);
% title('Satellite swath (optical)','FontSize',16,'interpreter','latex')
box on
saveLocation = 'D:\OneDrive - RMIT University\20. My Thesis\0. Images';
Filename = fullfile(saveLocation, 'Figure210');
% print(h_Fig, '-dpng','-r600',Filename)
%% Get ground reflectrivity - STEP3.Reflectivity Simulator Figure211
Scale = 0.9;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
pc =pcolor(xEast/1000,yNorth/1000,sigma);
% scatter(xEast(:)/1000,yNorth(:)/1000,2,sigma(:),'MarkerEdgeColor','none','MarkerFaceColor','flat')
ax=gca;
pc.LineStyle='none';
colormap bone
axis equal
hold on
line(xEast(1,:)/1000,yNorth(1,:)/1000,'LineStyle', '-', 'Color',ColorOrder(7,:), 'LineWidth', 2)
hold on
line(xEast(end,:)/1000,yNorth(end,:)/1000,'LineStyle', '-', 'Color',ColorOrder(7,:), 'LineWidth', 2)
hold on
Idx = round(length(xEast)/2);                                                                   % Index of mid point of the dwell
line(xEast(Idx,:)/1000,yNorth(Idx,:)/1000,'LineStyle', '--', 'Color',ColorOrder(5,:), 'LineWidth', 2)
hold on
plot(0,0,'+','LineWidth',1,'color',ColorOrder(7,:),'MarkerSize', 15);                           % Mid point (reference)
% xlabel('x-axis [km]')
% ylabel('y-axis [km]')
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',12);
% title('Satellite swath (optical)','FontSize',16,'interpreter','latex')
box on
% xticks([])
% yticks([])
xlim([-6, 6]);
grid on
% ylim([-4.5, 4.5]);
saveLocation = 'D:\OneDrive - RMIT University\20. My Thesis\0. Images';
Filename = fullfile(saveLocation, 'Figure220a');
print(h_Fig, '-dpng','-r600',Filename)
%% Plot Unfocused SAR raw data Figure5b
Scale = 1.1;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
pc =pcolor(FastTime/1e-6,1:etaTotal,abs(sqd));
pc.LineStyle='none';
colormap bone
% xlabel('Fast time [\mus]')
% ylabel('Azimuth index')
% title('Step 0: Raw time domain (magnitude)')
ax=gca;
pc.LineStyle='none';
grid on
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',12);
box on
% xticks([])
% yticks([])
saveLocation = 'D:\OneDrive - RMIT University\20. My Thesis\0. Images';
Filename = fullfile(saveLocation, 'Figure212');
% print(h_Fig, '-dpng','-r600',Filename)
%% Plot the stages of the compresion Figure212
Scale = 1.1;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
subplot(2,3,1)
pc =pcolor(FastTime/1e-6,1:etaTotal,real(sqd));xticks([]),yticks([])
pc.LineStyle='none';
colormap sky
xlabel('Fast time [\mus]')
ylabel('Azimuth index','FontSize',8)
title('Raw time domain','FontSize',8)
subplot(2,3,2)
plot(FastTime/1e-6,real(G));xticks([]),yticks([])
xlabel('Range Frequency [MHz]','FontSize',8)
ylabel('Real component','FontSize',8)
title('Range matched filter','FontSize',8)
subplot(2,3,3)
pc =pcolor(FastTime/1e-6,1:etaTotal,real(src));xticks([]),yticks([])
pc.LineStyle='none';
xlabel('Fast time [\mus]','FontSize',8)
ylabel('Azimuth index','FontSize',8)
title('Step 1: Range compression','FontSize',8)
subplot(2,3,4)
pc =pcolor(FastTime/1e-6,1:etaTotal,abs(S2));xticks([]),yticks([])
pc.LineStyle='none';
xlabel('Fast time [\mus]','FontSize',8)
ylabel('Azimuth index','FontSize',8)
title('Step 2: Azimuth FFT','FontSize',8)
subplot(2,3,5)
plot(1:etaTotal,DeltaR);xticks([]),yticks([])
xlabel('Azimuth index','FontSize',8)
ylabel('Range compensation [m]','FontSize',8)
title('Step 3.1: Range compensation profile','FontSize',8)
subplot(2,3,6)
pc =pcolor(FastTime/1e-6,1:etaTotal,real(S2));xticks([]),yticks([])
pc.LineStyle='none';
xlabel('Fast time [\mus]','FontSize',8)
ylabel('Azimuth index','FontSize',8)
title('Step 3.2: RCMC','FontSize',8)
ax=gca;
pc.LineStyle='none';
grid on
set(gca,'LooseInset',get(gca,'TightInset'));
box on
xticks([])
yticks([])
saveLocation = 'D:\OneDrive - RMIT University\20. My Thesis\0. Images';
Filename = fullfile(saveLocation, 'Figure212');
print(h_Fig, '-dpng','-r600',Filename)
%% Plot interfernce signal Figure215
Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
subplot(2,1,1)
spectrogram(signal_IQ,500,0,500,'','yaxis','centered')
title('QPSK modulated signal injected to the SAR signal');
subplot(2,1,2)
spectrogram(LoRaIQ,500,0,500,'','yaxis','centered')
title('LoRa signal injected to the SAR signal');
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',10);
saveLocation = 'D:\OneDrive - RMIT University\20. My Thesis\0. Images';
Filename = fullfile(saveLocation, 'Figure215');
% print(h_Fig, '-dpng','-r600',Filename1)
%% Plot LoRa Signal
% for plotting clear LoRa SF =9, BW = 0.8MHz
Scale = 0.9;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
spectrogram(LoRaIQ,500,0,500,RadPar.fs,'yaxis','centered')
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',10);
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',10);
% Filename1='Figure215';
% print(h_Fig, '-dpng','-r600',Filename1)
%% Plot Normalized PSD of SAR clean signal  Figure 302
co = colororder;
Scale = 0.7;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
% [SignalPSD , ff] = pwelch(sqdT(:),size(sqdIT,1),[],size(sqdIT,1),RadPar.fs,'centered');
[SignalPSD , ff] = pwelch(sqdT(:),size(sqdT,1),[],size(sqdT,1),RadPar.fs,'centered');
h = plot(ff/1e6, 10*log10(SignalPSD));
set(h, 'Color', 'k');
% pwelch(sqdIT(:),size(sqdIT,1),[],size(sqdIT,1),RadPar.fs,'centered')
xlabel('Frequency (MHz)');
ylabel('Power/Frequency (dB/Hz)');
ax=gca;hold on;
grid on;gx.Box = 'on';set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',12);
saveLocation = 'D:\OneDrive - RMIT University\20. My Thesis\0. Images';
Filename = fullfile(saveLocation, 'Figure5c');
print(h_Fig, '-dpng','-r600',Filename)
%% Plot Normalized PSD of SAR clean signal  Figure 1
co = colororder;
Scale = 0.7;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
% [SignalPSD , ff] = pwelch(sqdT(:),size(sqdIT,1),[],size(sqdIT,1),RadPar.fs,'centered');
[SignalPSD , ff] = pwelch(sqdT(:),size(sqdT,1),[],size(sqdT,1),RadPar.fs,'centered');
h = plot(ff/1e6, 10*log10(SignalPSD));
set(h, 'Color', 'k');
% pwelch(sqdIT(:),size(sqdIT,1),[],size(sqdIT,1),RadPar.fs,'centered')
xlabel('Frequency (MHz)');
ylabel('Power/Frequency (dB/Hz)');
% title('PSD of SAR raw data before RFI removal')
ax=gca;hold on;
RFIfreqshift = 7.9995e6;
% Draw vertical lines at peak locations
% xline(RFIfreqshift*1e-6,'--','color',co(7,:),'LineWidth',1);
% Show the SIR-PSD on the plot 
yline(PSD_sqd,'--','color', co(5,:),'LineWidth',1)
RFI_Pr_dB = 20*log10(rms(sInterfT(:)));                      
PSD_RFI = RFI_Pr_dB - 10*log10(IR.bw);            % Average PSD of RFI
% yline(PSD_RFI,'--','color', co(1,:),'LineWidth',1)
% Use the text function to place the text
x_text = -15.3; 
y_text = (PSD_sqd + PSD_RFI) /2;
set(gca, 'XTick', sort([RFIfreqshift*1e-6, get(gca, 'XTick')]));xtickformat('%,.0f');
ax.XTickLabel{7} = ['\color[rgb]{' num2str(139/255) ', ' num2str(0/255) ', ' num2str(0/255) '}' ax.XTickLabel{7}];
% xlim([-5 15])
legend('','Average SAR PSD','Location','southwest')
grid on;gx.Box = 'on';set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',12);
hold off
% Filename1='Figure5a';
saveLocation = 'D:\OneDrive - RMIT University\20. My Thesis\0. Images';
Filename = fullfile(saveLocation, 'Figure5c');
print(h_Fig, '-dpng','-r600',Filename)
