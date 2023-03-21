clc;clear;close all hidden;
% load('matlabTestRef')                               % This is for reconstructing the ground reference point 
% load('matlabTestTargets')                           % This is for reconstructing three targets test
% load('matlabOptical1')                              % This is for reconstructing the whole scene image
load('Test01.mat')                                  % This is for testing
%% This is a raw-wise FFT / IFFT
fft1d2 = @ (x) fftshift(fft(fftshift(x,2),[],2),2);
ifft1d2 = @ (x) ifftshift(ifft(ifftshift(x,2),[],2),2);

% This is a cloumn-wise FFT - Azimuth
fft1d1 = @ (x) fftshift(fft(fftshift(x,1),[],1),1);
ifft1d1 = @ (x) ifftshift(ifft(ifftshift(x,1),[],1),1);
%% STEP6.SAR Image Processing
%% Step 1: Range Compression
So = fft1d2(sqd);                                   % FFT the time domain signal (FFT along each eta row)

% This is the template for the Matched Filter
tau = 0; 
sb = exp(-1j*pi *   (2*RadPar.fo * tau - RadPar.K*(FastTime-tau).^2   )    ) ...
    .*(FastTime>(-RadPar.T/2+tau)).*(FastTime<(RadPar.T/2+tau));
G = conj(fftshift(fft(fftshift(sb))));

Src  = repmat(G,size(So,1),1).*So;                  % Equation 5.5 - Multiplying the filter by So (Frequency Domain Multiplication)
src  = ifft1d2(Src);                                % Equation 5.6 - Inverse Fourier transform along each pulse to put the data back into range
%% Generate the reference signal based on the ground referene point for the migration step
So_ref = fft1d2(sqd_ref);                           % FFT the time domain reference signal (FFT along each eta row)
Src_ref = repmat(G,size(So_ref,1),1).*So_ref;       % Equation 5.5 - Multiplying the filter by So (Frequency Domain Multiplication)
src_ref = ifft1d2(Src_ref);                         % Equation 5.6 - Inverse Fourier transform along each pulse to put the data back into range
%% Step 2 Azimuth FFT
S1_ref = fft1d1(src_ref);
S1 = fft1d1(src);
%% Step 3 Range cell migration compensation
DeltaR = R-Ro;

% Shifting range cells
RangeBin = 2*RadPar.ts*c;                           % Ground range resolution
NbinsShift = round(DeltaR/RangeBin);
for ctr=1:size(S1,1)
    S2(ctr,:) = circshift(S1(ctr,:),NbinsShift);
end
for ctr=1:size(S1_ref,1)
    S2_ref(ctr,:) = circshift(S1_ref(ctr,:),NbinsShift);
end
%% Step 4 Azimuth compression
%%% Analytical Method %%%%
% R = R - DeltaR;
% Haz = exp(1j*2*pi* R *2*RadPar.fo/c);             % Azimuth Matched Filter
% Haz = exp(-1j*2*pi* DeltaR *2*RadPar.fo/c);         % Azimuth Matched Filter
% S3 = S2 .* repmat(Haz,1,size(S2,2));              % Compressed data after azimuth compression

% FrequencyAzimuth = -Param.PRF/2 : Param.PRF/ size(sqd,1): Param.PRF/2-(Param.PRF/size(sqd,1));    % Slowtime frequency array - Azimuth frequency
% [~,velocity,~] = states(sat);
% v = abs(sqrt((velocity(1,:).^2 + velocity(2,:).^2+ velocity(3,:).^2)));
% Ka = 2 * v.^2 ./ (RadPar.Lambda *  DeltaR.');
% Haz = exp(1j*pi* FrequencyAzimuth.^2 ./Ka);         % Azimuth Matched Filter
% S3 = S2 .* repmat(Haz.',1,size(S2,2));              % Compressed data after azimuth compression

%%% Empirical Method %%%%
midpoint = round(size(S2_ref,2)/2)+1;
S2_ref = repmat(S2_ref(:,midpoint),1,size(S2,2));
Haz= conj(S2_ref);
S3 = S2 .* Haz;
%% Step 5 Azimuth IFFT
sSLC = ifft1d1(S3);                          % Final Focused SAR Image
%% Plot Focused SAR Image
figure(1)
Range =(-(numel(FastTime)/2)*RangeBin:RangeBin:(numel(FastTime)/2-1)*RangeBin);
speed= mean(sqrt(sum((diff(SatECI,[],2)).^2)) /Param.ts);
CrossRange = (1:etaTotal)*Param.ts*speed/1000;

sSLC=sSLC./max(abs(sSLC),[],"all");
ax=gca;
pc =pcolor(Range/1000,CrossRange,(abs(sSLC)));
pc.LineStyle='none';
ax.YAxis.Direction = 'reverse';
ax.XAxis.Direction = 'reverse';
xlabel('Sar Range [km]')
ylabel('Cross Range [km]')
title('Step 5: Comressed image')
colormap bone
axis equal
%% Geographic projection
% Here we map the cross-range / range to domain to the original geographic
% coordinates longitute / latitude.

for eta =1:etaTotal
[~,~,SAR_Range(eta,:)] = geodetic2aer(Targetlat(eta,:),Targetlon(eta,:),0,Satlla(eta,1),Satlla(eta,2),Satlla(eta,3),E);
SAR_Range(eta,:) = (SAR_Range(eta,:)-slantrangeMid(eta));
SARlat(eta,:) = interp1(SAR_Range(eta,:),Targetlat(eta,:),Range,"linear","extrap");
SARlon(eta,:) = interp1(SAR_Range(eta,:),Targetlon(eta,:),Range,"linear","extrap");
end
figure(2)
subplot(1,2,1)
ax=gca;
[xEast,yNorth,~] = latlon2local(SARlat,SARlon,0,GRP);
%scatter(xEast(:)/1000,yNorth(:)/1000,2,double(abs(sSLC(:))),'MarkerEdgeColor','none','MarkerFaceColor','flat')
J = imadjust(double(abs(sSLC)));
imagesc(xEast(1,:)/1000,yNorth(:,1)/1000,J)
colormap bone
ax.YAxis.Direction="Normal";
axis equal
hold on
plot(0,0,'+'); % Mid point (reference)
xlabel('x-axis [km]')
ylabel('y-axis [km]')
title('Processed SAR Image')

subplot(1,2,2)
ax=gca;
[xEast,yNorth,~] = latlon2local(Targetlat,Targetlon,0,GRP);
%scatter(xEast(:)/1000,yNorth(:)/1000,2,a(:),'MarkerEdgeColor','none','MarkerFaceColor','flat')
a=a./sum(a,"all");
J = imadjust(a);
imagesc(xEast(1,:)/1000,yNorth(:,1)/1000,J)
colormap bone
ax.YAxis.Direction="Normal";
axis equal
hold on
plot(0,0,'+'); % Mid point (reference)
xlabel('x-axis [km]')
ylabel('y-axis [km]')
title('Satellite swath (optical)')
