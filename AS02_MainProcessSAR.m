clc; clear; close all
load('Test01.mat')
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
S1_ref = fft1d1(src_ref);

%% Plotting the range-compressed image
subplot(2,4,2)
pc =pcolor(FastTime/1e-6,1:etaTotal,real(src));
pc.LineStyle='none';
xlabel('Fast time [\mus]')
ylabel('Azimuth index')
title('Step 1: Range compression')
drawnow
%% Step 2 Azimuth FFT
S1 = fft1d1(src);
subplot(2,4,3)
pc =pcolor(FastTime/1e-6,1:etaTotal,abs(S1));
pc.LineStyle='none';
xlabel('Fast time [ms]')
ylabel('Azimuth index')
title('Step 2: Az FFT')

%% Step 3 Range cell migration compensation
DeltaR = R - Ro;
subplot(2,4,4)
plot(1:etaTotal,abs(DeltaR));
xlabel('Azimuth index')
ylabel('Range compensation [m]')
title('Step 3.1: Range compensation profile')

% shiftring range cells
RangeBin = 2*RadPar.ts*c;
NbinsShift = round(DeltaR/RangeBin);
for ctr=1:size(S1,1)
    S2(ctr,:) = circshift(S1(ctr,:),NbinsShift);
end

for ctr=1:size(S1_ref,1)
    S2_ref(ctr,:) = circshift(S1_ref(ctr,:),NbinsShift);
end
subplot(2,4,5)
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
Range =(0:RangeBin:(numel(FastTime)-1)*RangeBin)/1000/2;
speed= mean(sqrt(sum((diff(SatECI,[],2)).^2)) /Param.dt);

CrossRange = (1:etaTotal)*Param.dt*speed/1000;
%%
sSLC=sSLC./max(abs(sSLC),[],"all");
ax=gca;
pc =pcolor(Range,CrossRange,(abs(sSLC)));
pc.LineStyle='none';
ax.YAxis.Direction = 'reverse';
ax.XAxis.Direction = 'reverse';
xlabel('Range [km]')
ylabel('Cross Range [km]')
title('Step 5: Comressed image')
colormap bone
axis equal

%%
%figure(3)
% subplot(2,1,1)
% clf
% plot(real(S2(:,round(length(S2)/2+2))))
% hold on
% plot(real(Haz)*1e-10)
% plot(imag(sSLC(:,round(length(sSLC)/2+2))))

