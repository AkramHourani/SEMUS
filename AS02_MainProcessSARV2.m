clc; clear; close all
close all hidden;
%load('Test01')
%load('Test02')
load('SAR_Image')
%% This is a raw-wise FFT / IFFT
fft1d2 = @ (x) fftshift(fft(fftshift(x,2),[],2),2);
ifft1d2 = @ (x) ifftshift(ifft(ifftshift(x,2),[],2),2);
% This is a cloumn-wise FFT - Azimuth
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
%% STEP6.SAR Image Processing
%% Step 1: Range Compression
So = fft1d2(sqd);                                   % FFT the time domain signal (FFT along each eta row)
% This is the template for the Match Filter
tau = 0;
sb = exp(1j*pi *   (-2*RadPar.fo * tau + RadPar.K*(FastTime-tau).^2   )    ) ...
    .*(FastTime>(-RadPar.T/2+tau)).*(FastTime<(RadPar.T/2+tau));
G = conj(fftshift(fft(fftshift(sb))));
Src  = repmat(G,size(So,1),1).*So;               % Equation 5.5 - Multiplying the filter by So (Frequency Domain Multiplication)
src  = ifft1d2(Src);                             % Equation 5.6 - Inverse Fourier transform along each pulse to put the data back into range
%% Generate the reference signal based on the ground referene point
So = fft1d2(sqd_ref);                            % FFT the time domain reference signal (FFT along each eta row)
% This is the template for the Match Filter
tau = 0;
sb = exp(1j*pi *   (-2*RadPar.fo * tau + RadPar.K*(FastTime-tau).^2   )    ) ...
    .*(FastTime>(-RadPar.T/2+tau)).*(FastTime<(RadPar.T/2+tau));
G = conj(fftshift(fft(fftshift(sb))));
Src_ref  = repmat(G,size(So,1),1).*So;          % Equation 5.5 - Multiplying the filter by So (Frequency Domain Multiplication)
src_ref  = ifft1d2(Src_ref);                    % Equation 5.6 - Inverse Fourier transform along each pulse to put the data back into range
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

% shifting the range cells
RangeBin = RadPar.ts*c/2;
NbinsShift = -round(DeltaR/RangeBin)*2;
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
Haz = exp(1j*pi*R*4*RadPar.fo/c);               % Azimuth Analytical Matched Filter
%S3 = S2 .* repmat(Haz,1,size(S2,2));
subplot(2,4,6)
plot(real(Haz))

%%% Empirical Method %%%%
midpoint = round(size(S2_ref,2)/2)+1;
S2_ref = repmat(S2_ref(:,midpoint),1,size(S2,2));
S3 = S2 .* conj(S2_ref);
%% Step 5 Azimuth IFFT
sSLC = ifft1d1(S3);                             % Final Focused SAR Image

%% plotting (this is an approzimate projection of the swath)
figure(2)
Range =(-(numel(FastTime)/2)*RangeBin:RangeBin:(numel(FastTime)/2-1)*RangeBin);
[~,El,~]= geodetic2aer(GRP(1),GRP(2),0,Satlla(1,1),Satlla(1,2),Satlla(1,3),E);

RangeGround = Range/cosd(El);

speed= mean(sqrt(sum((diff(SatECI,[],2)).^2)) /Param.ts);
CrossRange = (1:etaTotal)*Param.ts*speed/1000;

J = abs(sSLC);
J = J./max(J,[],"all");
J = imresize(J,2);
J  = imgaussfilt(J ,2); % Smothing filter
J = imadjust(J,[0 0.6]);
imagesc(RangeGround/1000,CrossRange,J)
ax=gca;
ax.YAxis.Direction = 'reverse';
ax.XAxis.Direction = 'reverse';
xlabel('Ground distance (approx.) [km]')
ylabel('Cross Range [km]')
title('Step 5: Comressed image')
colormap bone
axis equal
xlim([-1 1]*Swathwidth/2/1000);
