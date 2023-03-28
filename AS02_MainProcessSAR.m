clc; clear; close all
close all hidden;
load('SAR_Image.mat')
%% This is a raw-wise FFT / IFFT
fft1d2 = @ (x) fftshift(fft(fftshift(x,2),[],2),2);
ifft1d2 = @ (x) ifftshift(ifft(ifftshift(x,2),[],2),2);
% This is a cloumn-wise FFT - Azimuth
fft1d1 = @ (x) fftshift(fft(fftshift(x,1),[],1),1);
ifft1d1 = @ (x) ifftshift(ifft(ifftshift(x,1),[],1),1);
%% plotting raw time domain signal
figure(1);
subplot(2,3,1)
pc =pcolor(FastTime/1e-6,1:etaTotal,real(sqd));
pc.LineStyle='none';
colormap parula
xlabel('Fast time [\mus]')
ylabel('Azimuth index')
title('Step 0: Raw time domain (magnitude)')
%% STEP5.SAR Image Processing
%% Step 1: Range Compression
So = fft1d2(sqd);                                % FFT the time domain signal (FFT along each eta row)
% This is the template for the Match Filter
tau = 0;
sb = exp(1j*pi *   (-2*RadPar.fo * tau + RadPar.K*(FastTime-tau).^2   )    ) ...
    .*(FastTime>(-RadPar.T/2+tau)).*(FastTime<(RadPar.T/2+tau));
G = conj(fftshift(fft(fftshift(sb))));
Src  = repmat(G,size(So,1),1).*So;               % Equation 5.5 - Multiplying the filter by So (Frequency Domain Multiplication)
src  = ifft1d2(Src);                             % Equation 5.6 - Inverse Fourier transform along each pulse to put the data back into range
%% Generate the reference signal based on the ground referene point
So_ref = fft1d2(sqd_ref);                        % FFT the time domain signal (FFT along each eta row)
Src_ref  = repmat(G,size(So,1),1).*So_ref; 
src_ref  = ifft1d2(Src_ref);
%% Plotting the range-compressed image
subplot(2,3,2)
pc =pcolor(FastTime/1e-6,1:etaTotal,real(src));
pc.LineStyle='none';
xlabel('Fast time [\mus]')
ylabel('Azimuth index')
title('Step 1: Range compression')
drawnow
%% Step 2 Azimuth FFT
S2_ref = fft1d1(src_ref);
S2 = fft1d1(src);
subplot(2,3,3)
pc =pcolor(FastTime/1e-6,1:etaTotal,abs(S2));
pc.LineStyle='none';
xlabel('Fast time [ms]')
ylabel('Azimuth index')
title('Step 2: Azimuth FFT')
%% Step 3 Range cell migration compensation
DeltaR = R - Ro;                                    % Based on the GRP
subplot(2,3,4)
plot(1:etaTotal,DeltaR);
xlabel('Azimuth index')
ylabel('Range compensation [m]')
title('Step 3.1: Range compensation profile')

% Shiftring the range cells
RangeBin = RadPar.ts*c;
NbinsShift = -round(DeltaR/RangeBin);
for AzCtr=1:etaTotal
    S2(AzCtr,:) = circshift(S2(AzCtr,:),NbinsShift(AzCtr));
    S2_ref(AzCtr,:) = circshift(S2_ref(AzCtr,:),NbinsShift);
end

subplot(2,3,5)
pc =pcolor(FastTime/1e-6,1:etaTotal,real(S2));
pc.LineStyle='none';
xlabel('Fast time [\mus]')
ylabel('Azimuth index')
title('Step 3.2: RCMC')
drawnow
%% Step 4 Azimuth compression
Haz = exp(-1j*pi*R*4*RadPar.fo/c);                  % Azimuth Analytical Matched Filter
%S3 = S2 .* repmat(Haz,1,size(S2,2));
subplot(2,3,6)
plot(real(Haz))

%%% Empirical Method %%%%
midpoint = round(size(S2_ref,2)/2)+1;
S2_ref = repmat(S2_ref(:,midpoint),1,size(S2,2));
S3 = S2 .* conj(S2_ref);
%% Step 5 Azimuth IFFT
sSLC = ifft1d1(S3);                                 % Final Focused SAR Image
%% Plotting (An approximate projection of the swath)
figure(2)
clf
%speed= mean(sqrt(sum((diff(SatECI,[],2)).^2)) /Param.ts);
%CrossRange = (1:etaTotal)*Param.ts*speed/1000;
Img=abs(sSLC)./max(abs(sSLC),[],"all");
Calibration = 1;

% Time equivalent range (i.e. twice the slant range in case of mono-staitic SAR
RangeEq =(-(numel(FastTime)/2+Calibration)*RangeBin:RangeBin:(numel(FastTime)/2-Calibration-1)*RangeBin);
ax=gca;
pc =pcolor(RangeEq/1000,1:size(Img,1),Img);
pc.LineStyle='none';
ax.YAxis.Direction = 'reverse';
ax.XAxis.Direction = 'reverse';
xlabel('Slant Range [km]')
ylabel('Cross Range [km]')
title('Step 5: Compressed image')
colormap turbo
%axis equal
drawnow
%imwrite(abs(sSLC),"Mono_SAR.png");
%% Geographic projection
%% First: Create transformation control points in Lat/Lon domain
ResAz = 10;
ResR  = 10;
if etaTotal> ResAz
    n = round(etaTotal/ResAz);
    etaVec = downsample(1:etaTotal,n);              % Take a subsample from the slow time vector
    ResAz = length(etaVec);
else 
    ResAz=etaTotal;
end

for AzCtr =1:ResAz % Creat the points based on the swath
    eta =etaVec(AzCtr);
    % build the conneting line between the swath points
    [CLat(AzCtr,:),CLon(AzCtr,:)] = gcwaypts(latSwathL1(eta),lonSwathL1(eta),latSwathL2(eta),lonSwathL2(eta),ResR-1);

end
%% Second: calculate the transformation control points in Az/Range domain based on given Lat/Lon
for Ctr=1:length(CLat(:))
    % Find the disance profile of the target over the entire aquisition period
    [~,~,RTemp]=geodetic2aer(CLat(Ctr),CLon(Ctr),0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);
    [row,col] = ind2sub(size(CLat),Ctr);        % Find the row column indices of the control point
    % Find the closest approach    
    [Rmin, etaMin] = min(RTemp); 
    ContAz(row,col) = etaMin;                   % This is the slow-time of the closest approach
    ContR(row,col)  = (Rmin-Ro)*2;              % This is the SAR distance to the control point
end

% Third: Create the inverse transmofation function Az/Range -> Lat/Lon
AzR2LatLon = fitgeotform2d([CLat(:),CLon(:)],[ContAz(:),ContR(:)],"polynomial",4);
LatLon2AzR = fitgeotform2d([ContAz(:),ContR(:)],[CLat(:),CLon(:)],"polynomial",4);

% Test the tranformation
Output1 = transformPointsInverse(LatLon2AzR,[Targetlat(:), Targetlon(:)]);
hold on
%scatter(Output1(:,2)/1000,Output(:,1),"+","MarkerEdgeColor",ax.ColorOrder(2,:))
%% Third: Transform to local cartisian coordiantes
figure(3)
clf
% Transform from AzR -> Lat/Lon
etaVecM = repmat((1:etaTotal)',1,length(RangeEq));
RangeM  = repmat(RangeEq,etaTotal,1);
[LatImg, LonImg] = transformPointsInverse(AzR2LatLon,etaVecM, RangeM);

% Transfrom from Lat/Lon to NEC
[xImg,yImg,~]     = latlon2local(LatImg,LonImg,0,GRP);
%scatter(xImg(:),yImg(:),"+","MarkerEdgeColor",ax.ColorOrder(2,:))
hold on
axis equal
pc =pcolor(xImg/1000,yImg/1000,Img);
%scatter(xEast(:)/1000,yNorth(:)/1000,"o","MarkerEdgeColor",ax.ColorOrder(2,:))
scatter(0,0,"+","MarkerEdgeColor",ax.ColorOrder(1,:))
pc.LineStyle='none';
grid on
% colormap turbo
colormap bone
ax.YAxis.Direction = 'reverse';
ax.XAxis.Direction = 'reverse';
xlabel('North-axis [km]')
ylabel('East axis [km]')
title('Corrected geo image')
