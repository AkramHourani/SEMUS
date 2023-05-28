clc; clear; close all
%load('Mesh_Bi.mat')
load('Sydney.mat')
close all hidden;
ttlsz = 8;

%% This is a raw-wise FFT / IFFT

% fft performs a DFT using a fast (FFT) algorithm
% fftshift swaps the left and right portions of a sequence, placing the first element
% in the center
% ifft performs a IDFT using an inverse fast (FFT) algorithm
% ifftshift swaps the left and right portions of a sequence, placing the central
% element at the beginning

fft1d2 = @ (x) fftshift(fft(fftshift(x,2),[],2),2);
ifft1d2 = @ (x) ifftshift(ifft(ifftshift(x,2),[],2),2);
fft1d1 = @ (x) fftshift(fft(fftshift(x,1),[],1),1);
ifft1d1 = @ (x) ifftshift(ifft(ifftshift(x,1),[],1),1);

Scale = 1;
%h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[1000 150 800 800/1.618*Scale]);
h_Fig = figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 2*3.5*2 3.5*2/1.618*Scale],'Position',[1000 150 800 800/1.618*Scale]);

t = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact');
nexttile
%% plotting raw time domain signal
%subplot(2,3,1)
%pc = pcolor(FastTime/1e-6,1:etaTotal,real(sqd));
%pc.LineStyle='none';
imagesc(real(sqd))
colormap bone
xlabel('Fast time [\mus]')
ylabel('Azimuth index')
%xticklabels('')
%yticklabels('')
%xticks([])
%yticks([])

title('(a)', 'Units', 'normalized', 'Position',  [0.05,0.4,0],'interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','none');
%subtitle('Raw time domain (magnitude)')
set(gca,'LooseInset',get(gca,'TightInset'));

%% Step 1: Range Compression
So = fft1d2(sqd); % FFT the time domain signal (FFT along each eta row)

% This is the template for the Match Filter
tau = 0;
sb = exp(1j*pi *   (-2*RadPar.fo * tau + RadPar.K*(FastTime-tau).^2   )    ) ...
    .*(FastTime>(-RadPar.T/2+tau)).*(FastTime<(RadPar.T/2+tau));
G = conj(fftshift(fft(fftshift(sb))));
Src  = repmat(G,size(So,1),1).*So; %Equation 5.5
src  = ifft1d2(Src); %Equation 5.6

%% Plotting the Real component of the range matched filter
%subplot(2,3,2)
nexttile
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

%% Generate the reference signal based on the ground reference point
So = fft1d2(sqd_ref); % FFT the time domain signal (FFT along each eta row)
% This is the template for the Match Filter : complex conjugate of the
% input signal
tau = 0;
sb = exp(1j*pi *   (-2*RadPar.fo * tau + RadPar.K*(FastTime-tau).^2   )    ) ...
    .*(FastTime>(-RadPar.T/2+tau)).*(FastTime<(RadPar.T/2+tau));
G = conj(fftshift(fft(fftshift(sb))));  % G(f)
Src  = repmat(G,size(So,1),1).*So; %Equation 5.5
src_ref  = ifft1d2(Src); %Equation 5.6 

%% Plotting the range-compressed image
%subplot(2,3,3)
nexttile
%pc = pcolor(FastTime/1e-6,1:etaTotal,real(src));
imagesc(real(src));
%pc.LineStyle='none';
xlabel('Fast time [\mus]')
ylabel('Azimuth index')
%xticklabels('')
%yticklabels('')
%xticks([])
%yticks([])
title('(c)', 'Units', 'normalized', 'Position',  [0.05,0.4,0],'interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','none');
%subtitle('Step 1: Range compression')
set(gca,'LooseInset',get(gca,'TightInset'));

%drawnow

%% Step 2 Azimuth FFT
%subplot(2,3,4)
nexttile
S2_ref = fft1d1(src_ref);
S2 = fft1d1(src);
%subplot(2,4,5)
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

%% Step 3 Range cell migration compensation
%subplot(2,3,5)
nexttile
DeltaR = RSoI+RI-Ro; % The migration of the GRP
%subplot(2,4,6)
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

% shiftring range cells (Needs to be update for Bi Static SAR)
RangeBin = RadPar.ts*c/2;
NbinsShift = -round(DeltaR/RangeBin)*2;
% RangeBin = RadPar.ts*c;
% NbinsShift = -round(DeltaR/RangeBin);
for eta=1:etaTotal
    S2(eta,:) = circshift(S2(eta,:),NbinsShift(eta));
    S2_ref(eta,:) = circshift(S2_ref(eta,:),NbinsShift(eta));
end

%subplot(2,3,6)
nexttile
%pc = pcolor(FastTime/1e-6,1:etaTotal,real(src));
imagesc(real(src))
%pc.LineStyle='none';
xlabel('Fast time [\mus]')
ylabel('Azimuth index')
%xticklabels('')
%yticklabels('')
%xticks([])
%yticks([])

title('(f)', 'Units', 'normalized', 'Position',  [0.05,0.4,0],'interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','none');
%subtitle('Step 3.2: RCMC')
set(gca,'LooseInset',get(gca,'TightInset'));

drawnow
colormap bone
FilenameP5='Figure7';
print(h_Fig, '-dpng','-r600',FilenameP5)
%movefile([FilenameP5 '.png'],'Figures')

%% Step 4/5 Azimuth compression

%%% Empirical Method %%%%
midpoint = round(size(S2_ref,2)/2)+1;
S2_ref = repmat(S2_ref(:,midpoint),1,size(S2,2));
S3 = S2 .* conj(S2_ref);

%% Step 5 Azimuth IFFT
sSLC = ifft1d1(S3);

%% Plotting Compressed Image
%CalibrationR = 2;
%CalibrationAz = 160;
%Img = abs(sSLC)./max(abs(sSLC),[],"all");

Scale = 1;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[1000 150 800 800/1.618*Scale]);

Range =(-(numel(FastTime)/2)*RangeBin:RangeBin:(numel(FastTime)/2-1)*RangeBin);
[~,El,~]= geodetic2aer(GRP(1),GRP(2),0,SatllaSoI(1,1),SatllaSoI(1,2),SatllaSoI(1,3),E); %% I changed Satlla to SatllaSoI

RangeGround = Range/cosd(El);

speed = mean(sqrt(sum((diff(SatECISoI,[],2)).^2)) /Param.ts);%% I changed SatECI to SatECISoI
CrossRange = (1:etaTotal)*Param.ts*speed/1000;

J = abs(sSLC);
J = J./max(J,[],"all");
J = imresize(J,2);
J = imgaussfilt(J ,2); % Smoothing filter
J = imadjust(J,[0 0.6]);

imagesc(RangeGround/1000,CrossRange,J)

ax=gca;
ax.YAxis.Direction = 'reverse';
ax.XAxis.Direction = 'reverse';
xlabel('Ground distance (approx.) [km]')
ylabel('Azimuth [km]')
title('Compressed image')
colormap bone
axis equal
xlim([-1 1]*SwathwidthSoI/2/1000);%% I changed Swathwidth to SwathwidthSoI

% ax = gca;
% %pc = pcolor(Range/1000,etaVec,Img);
% %pc.LineStyle='none';
% ax.YAxis.Direction = 'reverse';
% ax.XAxis.Direction = 'reverse';
% % xticklabels('')
% % yticklabels('')
% xticks('')
% yticks('')
% xlabel('Fast-time','interpreter','Tex')
% ylabel('Azimuth Index','interpreter','Tex')
% title('Final : Compressed Image','interpreter','latex') % Compressed SAR Data, Single Look Complex
% set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',12);
% %colormap turbo
% colormap bone
% % axis equal
% drawnow
%
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',12);
FilenameP6='Figure9b';
print(h_Fig, '-dpng','-r600',FilenameP6)
%movefile([FilenameP6 '.png'],'Figures')

%%

% %% Geographic projection
% % First: Create transformation control points in Lat/Lon domain
% ResAz = 8;
% ResR  = 8;
% if etaTotal> ResAz
%     n = round(etaTotal/ResAz);
%     etaVec = downsample(1:etaTotal,n); % Take a subsample from the slow time vector
%     ResAz = length(etaVec);
% else 
%     ResAz=etaTotal;
% end
% CLat=[];
% CLon=[];
% disp("Creating control points in Lat/Lon...")
% for AzCtr =1:ResAz % Create the points based on the swath
%     eta = etaVec(AzCtr);
%     % build the connecting line between the swath points
%     [CLat(AzCtr,:),CLon(AzCtr,:)] = gcwaypts(latSwathL1SoI(eta),lonSwathL1SoI(eta),latSwathL2SoI(eta),lonSwathL2SoI(eta),ResR-1);
% end
% 
% % Second: claculate the transformation control points in Az/Range domain based on given Lat/Lon
% disp("Calculating transfomation points Lat/Lon -> AzRange...")
% ContAz=[];
% ContR=[];
% for Ctr=1:length(CLat(:))
%     % Find the distance profile of the target over the entire aquisition period
%     [~,~,RTemp1]=geodetic2aer(CLat(Ctr),CLon(Ctr),0,SatllaSoI(:,1),SatllaSoI(:,2),SatllaSoI(:,3),E);
%     [~,~,RTemp2]=geodetic2aer(CLat(Ctr),CLon(Ctr),0,SatllaI(:,1),SatllaI(:,2),SatllaI(:,3),E);
%     [row,col] = ind2sub(size(CLat),Ctr);
%    
%     % find the closest approach
%     [Rmin, etaMin]  = min(RTemp1+RTemp2); 
%     ContAz(row,col) = etaMin; % This is the slow-time of the closest approach
%     ContR(row,col)  = Rmin-Ro; % This is the slant distance to the control point
% end
% 
% % Third: Create the inverse transformation function Az/Range -> Lat/Lon
% AzR2LatLon = fitgeotform2d([CLat(:),CLon(:)],[ContAz(:),ContR(:)],"polynomial",4);
% LatLon2AzR = fitgeotform2d([ContAz(:),ContR(:)],[CLat(:),CLon(:)],"polynomial",4);
% 
% % Test the tranformation
% Output1 = transformPointsInverse(LatLon2AzR,[Targetlat(:), Targetlon(:)]);
% hold on
% 
% %%To compare initial targets with compressed image
% %scatter(Output1(:,2)/1000,Output1(:,1),"+","MarkerEdgeColor",ax.ColorOrder(2,:))
% 
% 
% %scatter(Output2(:,2),Output2(:,1))
% %scatter(CLon(:),CLat(:))
% 
% 
% 
% %%
% % Transform to local cartesian coordinates
% Scale = 1.1;
% h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[1000 150 800 800/1.618*Scale]);
% 
% clf
% % Transform from AzR -> Lat/Lon
% etaVec = CalibrationAz:1:(size(Img,1)+CalibrationAz-1);
% etaVecM = repmat(etaVec',1,length(Range));
% RangeM  = repmat(Range,etaTotal,1);
% [LatImg, LonImg] = transformPointsInverse(AzR2LatLon,etaVecM, RangeM);
% 
% % Transfrom from Lat/Lon to NEC
% [xImg,yImg,~]     = latlon2local(LatImg,LonImg,0,GRP);
% %scatter(xImg(:)/1000,yImg(:)/1000,"+","MarkerEdgeColor",ax.ColorOrder(2,:))
% hold on
% 
% pc = pcolor(xImg/1000,yImg/1000,Img);
% %scatter(xEast(:)/1000,yNorth(:)/1000,"+","MarkerEdgeColor",ax.ColorOrder(2,:)) 
% 
% pc.LineStyle='none';
% colormap turbo
% ax.YAxis.Direction = 'reverse';
% ax.XAxis.Direction = 'reverse';
% %xlim([-1 1]/2* SwathwidthSoI/1000)
% axis equal
% xlabel('x-axis [km]','interpreter','latex')
% ylabel('y-axis [km]','interpreter','latex')
% title('Final Step: Projected Image','interpreter','latex')
% set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',12);
% FilenameP7='FigureP7';
% print(h_Fig, '-dpng','-r600',FilenameP7)
