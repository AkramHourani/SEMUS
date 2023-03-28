clc; clear; close all
load('Mesh_Bi.mat')
%load('Sydney_Bi.mat')
close all hidden;
%% This is a raw-wise FFT / IFFT
fft1d2 = @ (x) fftshift(fft(fftshift(x,2),[],2),2);
ifft1d2 = @ (x) ifftshift(ifft(ifftshift(x,2),[],2),2);
fft1d1 = @ (x) fftshift(fft(fftshift(x,1),[],1),1);
ifft1d1 = @ (x) ifftshift(ifft(ifftshift(x,1),[],1),1);
%% plotting raw time domain signal
subplot(2,4,1)
pc =pcolor(FastTime/1e-6,1:etaTotal,real(sqd));
pc.LineStyle='none';
colormap parula
xlabel('Fast time [\mus]')
ylabel('Azimuth index')
title('Step 0: Raw time domain (magnitude)')
%% Step 1: Range Compression
So = fft1d2(sqd); % FFT the time domain signal (FFT along each eta row)
% This is the template for the Match Filter
tau = 0;
sb = exp(1j*pi *   (-2*RadPar.fo * tau + RadPar.K*(FastTime-tau).^2   )    ) ...
    .*(FastTime>(-RadPar.T/2+tau)).*(FastTime<(RadPar.T/2+tau));
G = conj(fftshift(fft(fftshift(sb))));
Src  = repmat(G,size(So,1),1).*So; %Equation 5.5
src  = ifft1d2(Src); %Equation 5.6
%% Generate the reference signal based on the ground referene point
So = fft1d2(sqd_ref); % FFT the time domain signal (FFT along each eta row)
% This is the template for the Match Filter
tau = 0;
sb = exp(1j*pi *   (-2*RadPar.fo * tau + RadPar.K*(FastTime-tau).^2   )    ) ...
    .*(FastTime>(-RadPar.T/2+tau)).*(FastTime<(RadPar.T/2+tau));
G = conj(fftshift(fft(fftshift(sb))));
Src  = repmat(G,size(So,1),1).*So; %Equation 5.5
src_ref  = ifft1d2(Src); %Equation 5.6

%% Plotting the range-compressed image
subplot(2,4,2)
pc =pcolor(FastTime/1e-6,1:etaTotal,real(src));
pc.LineStyle='none';
xlabel('Fast time [\mus]')
ylabel('Azimuth index')
title('Step 1: Range compression')
drawnow
%% Step 2 Azimuth FFT
S2_ref = fft1d1(src_ref);
S2 = fft1d1(src);
subplot(2,4,5)
pc =pcolor(FastTime/1e-6,1:etaTotal,abs(S2));
pc.LineStyle='none';
xlabel('Fast time [ms]')
ylabel('Azimuth index')
title('Step 2: Az FFT')

%% Step 3 Range cell migration compensation

DeltaR = RSoI+RI - (Ro); % The migration of the GRP
subplot(2,4,4)
plot(1:etaTotal,abs(DeltaR));
xlabel('Azimuth index')
ylabel('Range compensation [m]')
title('Step 3.1: Range compensation profile')

% shiftring range cells
RangeBin = RadPar.ts*c;
NbinsShift = -round(DeltaR/RangeBin);
for eta=1:etaTotal
    S2(eta,:) = circshift(S2(eta,:),NbinsShift(eta));
    S2_ref(eta,:) = circshift(S2_ref(eta,:),NbinsShift);
end

subplot(2,4,3)
pc =pcolor(FastTime/1e-6,1:etaTotal,real(src));
pc.LineStyle='none';
xlabel('Fast time [\mus]')
ylabel('Azimuth index')
title('Step 3.2: RCMC')
drawnow



%% Step 4/5 Azimuth compression
Haz = exp(1j*pi*(RSoI+RI)*4*RadPar.fo/c);
subplot(2,4,6)
plot(real(Haz))

%S3 = S2 .* repmat(Haz,1,size(S2,2));
midpoint = round(size(S2_ref,2)/2)+1;
S2_ref = repmat(S2_ref(:,midpoint),1,size(S2,2));
S3 = S2 .* conj(S2_ref);
sSLC = ifft1d1(S3);


figure(2)

%%
figure(2)
clf
speed= mean(sqrt(sum((diff(SatECISoI,[],2)).^2)) /Param.ts);
CrossRange = (1:etaTotal)*Param.ts*speed/1000;
Img=abs(sSLC)./max(abs(sSLC),[],"all");
CalibrationR = 2;
CalibrationAz = 160;
Range =(-(numel(FastTime)/2+CalibrationR)*RangeBin:RangeBin:(numel(FastTime)/2-CalibrationR-1)*RangeBin);
etaVec = CalibrationAz:1:(size(Img,1)+CalibrationAz-1);
%etaVecM = repmat(etaVec',1,length(Range));
%RangeM  = repmat(Range,etaTotal,1);
ax=gca;
pc =pcolor(Range/1000,etaVec,Img);
pc.LineStyle='none';
ax.YAxis.Direction = 'reverse';
ax.XAxis.Direction = 'reverse';
xlabel('Bistatic Range [km]')
ylabel('Az Index')
title('Step 5: Comressed image')
colormap turbo
% axis equal
drawnow
%
% Geographic projection

% First: Creat transformation control points in Lat/Lon domain
ResAz = 10;
ResR  = 10;
if etaTotal> ResAz
    n = round(etaTotal/ResAz);
    etaVec = downsample(1:etaTotal,n); % Take a subsample from the slow time vector
    ResAz = length(etaVec);
else 
    ResAz=etaTotal;
end
CLat=[];
CLon=[];
disp("Creating control points in Lat/Lon...")
for AzCtr =1:ResAz % Creat the points based on the swath
    eta =etaVec(AzCtr);
    % build the conneting line between the swath points
    [CLat(AzCtr,:),CLon(AzCtr,:)] = gcwaypts(latSwathL1SoI(eta),lonSwathL1SoI(eta),latSwathL2SoI(eta),lonSwathL2SoI(eta),ResR-1);

end

% Second: claculate the transformation control points in Az/Range domain based on given Lat/Lon
disp("Calulating transfomation points Lat/Lon -> AzRange...")
ContAz=[];
ContR=[];
for Ctr=1:length(CLat(:))
    % Find the disance profile of the target over the entire aquisition period
    [~,~,RTemp1]=geodetic2aer(CLat(Ctr),CLon(Ctr),0,SatllaSoI(:,1),SatllaSoI(:,2),SatllaSoI(:,3),E);
    [~,~,RTemp2]=geodetic2aer(CLat(Ctr),CLon(Ctr),0,SatllaI(:,1),SatllaI(:,2),SatllaI(:,3),E);
    [row,col] = ind2sub(size(CLat),Ctr);
   
    % find the closest approach
    [Rmin, etaMin]  = min(RTemp1+RTemp2); 
    ContAz(row,col) = etaMin; % This is the slow-time of the closest approach
    ContR(row,col)  = Rmin-Ro; % This is the slant distance to the control point
    
end

% Third: Create the inverse transmofation function Az/Range -> Lat/Lon
AzR2LatLon = fitgeotform2d([CLat(:),CLon(:)],[ContAz(:),ContR(:)],"polynomial",4);
LatLon2AzR = fitgeotform2d([ContAz(:),ContR(:)],[CLat(:),CLon(:)],"polynomial",4);

% Test the tranformation
Output1 = transformPointsInverse(LatLon2AzR,[Targetlat(:), Targetlon(:)]);
hold on
scatter(Output1(:,2)/1000,Output1(:,1),"+","MarkerEdgeColor",ax.ColorOrder(2,:))
%scatter(Output2(:,2),Output2(:,1))
%scatter(CLon(:),CLat(:))
%%
% Transform to local cartisian coordiantes
figure(3)
clf
% Transform from AzR -> Lat/Lon
etaVec = CalibrationAz:1:(size(Img,1)+CalibrationAz-1);
etaVecM = repmat(etaVec',1,length(Range));
RangeM  = repmat(Range,etaTotal,1);
[LatImg, LonImg] = transformPointsInverse(AzR2LatLon,etaVecM, RangeM);

% Transfrom from Lat/Lon to NEC
[xImg,yImg,~]     = latlon2local(LatImg,LonImg,0,GRP);
%scatter(xImg(:)/1000,yImg(:)/1000,"+","MarkerEdgeColor",ax.ColorOrder(2,:))
hold on

pc =pcolor(xImg/1000,yImg/1000,Img);
%scatter(xEast(:)/1000,yNorth(:)/1000,"+","MarkerEdgeColor",ax.ColorOrder(2,:))

pc.LineStyle='none';
colormap turbo
ax.YAxis.Direction = 'reverse';
ax.XAxis.Direction = 'reverse';
%xlim([-1 1]/2* SwathwidthSoI/1000)
axis equal
xlabel('x-axis [km]')
ylabel('y-axis [km]')

