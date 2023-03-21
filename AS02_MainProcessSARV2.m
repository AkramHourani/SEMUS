clc; clear; close all
load('Test01.mat')
%load('matlabOptical.mat')
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
%% Step 3 Range cell migration compensation
DeltaR = R - Ro; % Based on the GRP
subplot(2,4,4)
plot(1:etaTotal,abs(DeltaR));
xlabel('Azimuth index')
ylabel('Range compensation [m]')
title('Step 3.1: Range compensation profile')

% shiftring range cells
RangeBin = RadPar.ts*c/2;
NbinsShift = -round(DeltaR/RangeBin);
for ctr=1:size(src,1)
    src(ctr,:) = circshift(src(ctr,:),NbinsShift(ctr));
end

for ctr=1:size(src_ref,1)
    src_ref(ctr,:) = circshift(src_ref(ctr,:),NbinsShift(ctr));
end
subplot(2,4,3)
pc =pcolor(FastTime/1e-6,1:etaTotal,real(src));
pc.LineStyle='none';
xlabel('Fast time [\mus]')
ylabel('Azimuth index')
title('Step 3.2: RCMC')
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
Param.ts = Param.dt

%% Step 4/5 Azimuth compression
Haz = exp(1j*pi*R*4*RadPar.fo/c);
subplot(2,4,6)
plot(real(Haz))

%S3 = S2 .* repmat(Haz,1,size(S2,2));
midpoint = round(size(S2_ref,2)/2)+1;
S2_ref = repmat(S2_ref(:,midpoint),1,size(S2,2));
S3 = S2 .* conj(S2_ref);
sSLC = ifft1d1(S3);


figure(2)
Range =(-(numel(FastTime)/2)*RangeBin:RangeBin:(numel(FastTime)/2-1)*RangeBin);
speed= mean(sqrt(sum((diff(SatECI,[],2)).^2)) /Param.ts);

CrossRange = (1:etaTotal)*Param.ts*speed/1000;
%%
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

if length(Targetlat)>5
    for eta =1:etaTotal
        [~,~,SAR_Range(eta,:)] = geodetic2aer(Targetlat(eta,:),Targetlon(eta,:),0,Satlla(eta,1),Satlla(eta,2),Satlla(eta,3),E);
        SAR_Range(eta,:) = (SAR_Range(eta,:)-slantrangeMid(eta));
        RangeExact = linspace(min(SAR_Range(eta,:)),max(SAR_Range(eta,:)),numel(FastTime));
        SARlat(eta,:) = interp1(SAR_Range(eta,:),Targetlat(eta,:),RangeExact,"linear","extrap");
        SARlon(eta,:) = interp1(SAR_Range(eta,:),Targetlon(eta,:),RangeExact,"linear","extrap");
    end
    figure
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
end