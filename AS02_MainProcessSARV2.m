clc; clear; close all
load('Mesh_Mono.mat')
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
DeltaR = R- Ro; % Based on the GRP
% shifting range cells (RCMC)
% RangeBin = RadPar.ts*c/2;
% NbinsShift = -round(DeltaR/RangeBin);
% subplot(2,4,3)
% plot(1:etaTotal,NbinsShift);
% xlabel('Azimuth index')
% ylabel('Range compensation [m]')
% title('Step 3.1: Range compensation profile')
% 
% 
% for eta=1:etaTotal
% src(eta,:) = circshift(src(eta,:),NbinsShift(eta));
% src_ref(eta,:) = circshift(src_ref(eta,:),NbinsShift);
% end

% subplot(2,4,4)
% pc =pcolor(FastTime/1e-6,1:etaTotal,real(src));
% pc.LineStyle='none';
% xlabel('Fast time [\mus]')
% ylabel('Azimuth index')
% title('RCMC')
% drawnow

%% Step 2 Azimuth FFT
S2_ref = fft1d1(src_ref);
S2 = fft1d1(src);
subplot(2,4,6)
pc =pcolor(FastTime/1e-6,1:etaTotal,abs(S2));
pc.LineStyle='none';
xlabel('Fast time [ms]')
ylabel('Azimuth index')
title('Step 2: Az FFT')


%% Step 3 Range cell migration corretion
DeltaR = R- Ro; % Based on the GRP
subplot(2,4,4)
plot(1:etaTotal,DeltaR);
xlabel('Azimuth index')
ylabel('Range compensation [m]')
title('Step 3.1: Range compensation profile')

% shiftring range cells
RangeBin = RadPar.ts*c/2;
NbinsShift = -floor(DeltaR/RangeBin);
for eta=1:etaTotal
S2(eta,:) = circshift(S2(eta,:),NbinsShift(eta));
S2_ref(eta,:) = circshift(S2_ref(eta,:),NbinsShift);
end


subplot(2,4,3)
pc =pcolor(FastTime/1e-6,1:etaTotal,real(S2));
pc.LineStyle='none';
xlabel('Fast time [\mus]')
ylabel('Azimuth index')
title('Step 3.2: RCMC')
drawnow

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
xlabel('Slant Range [km]')
ylabel('Cross Range [km]')
title('Step 5: Comressed image')
colormap bone
axis equal

%% Geographic projection
% Here we map the cross-range / range to domain to the original geographic
% coordinates longitute / latitude.

if length(Targetlat)>100
    for Ctr =1:length(Targetlat)
        [~,~,SAR_Range(Ctr,:)] = geodetic2aer(Targetlat(Ctr,:),Targetlon(Ctr,:),0,Satlla(Ctr,1),Satlla(Ctr,2),Satlla(Ctr,3),E);
        SAR_Range(Ctr,:) = (SAR_Range(Ctr,:)-slantrangeMid(Ctr));
       
        % Resample the range vector to match the number of points in the
        % fast time vector
        x = 1:length(SAR_Range(Ctr,:));
        xq = linspace(1,length(SAR_Range(Ctr,:)),length(FastTime));
        RangeExact = interp1(x,SAR_Range(Ctr,:),xq);;

        SARlat(Ctr,:) = interp1(SAR_Range(Ctr,:),Targetlat(Ctr,:),RangeExact,"linear","extrap");
        SARlon(Ctr,:) = interp1(SAR_Range(Ctr,:),Targetlon(Ctr,:),RangeExact,"linear","extrap");
    end
    figure
    subplot(1,2,1)
    ax=gca;
    [xEast,yNorth,~] = latlon2local(SARlat,SARlon,0,GRP);
    %scatter(xEast(:)/1000,yNorth(:)/1000,2,double(abs(sSLC(:))),'MarkerEdgeColor','none','MarkerFaceColor','flat')
    J = double(abs(sSLC));
    J = J/max(J,[],"all");
    J  = imgaussfilt(J ,2); % Smothing filter
    J = imadjust(J);
    
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
    a=a./max(a,[],"all");
    %J = imadjust(a);
    imagesc(xEast(1,:)/1000,yNorth(:,1)/1000,a)
    colormap bone
    ax.YAxis.Direction="Normal";
    axis equal
    hold on
    plot(0,0,'+'); % Mid point (reference)
    xlabel('x-axis [km]')
    ylabel('y-axis [km]')
    title('Satellite swath (optical)')
end