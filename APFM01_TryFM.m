clc
clear 
close all hidden

%% 0. Load parameters
%AP00_ParametersFerdi
Param.NtargetsAz = 11; % number of targets in each eta bin
Param.NtargetsRange = 11; % number of targets in each eta bin

% Frequency modulation : S(t)=Accos(2πfct+β∫m(t)dt)
% Pr = (Pt*(G^2)*(λ^2)*σ)/(((4π)^3)*(R^4))
% (SNR)_o = (Pt*(G^2)*(λ^2)*σ)/(((4π)^3)*k*Ts*B*L*(R^4))

%% FM Transmitter
fc      = 100e6; % Carrier freq = 100 MHz
fm      = 1;
beta    = 5;
tp      = 0:0.001:10;
message = cos(2*pi*fm*tp); %modulating signal, fm is message frequency
fm      = cos(2*pi*fc*tp+beta.*sin(2*pi*fm*tp)); % FM is modulated signal

figure(1)
subplot(3,1,1)
plot(tp, message);
title('Original signal');
xlabel('Time');
ylabel('Amplitude')

subplot(3,1,2)
plot(tp, fm);
title('FM signal');
xlabel('Time');
ylabel('Amplitude')

%% Radar Parameters
pt          = 100; % peak power in Watts
freq        = fc; % radar operating frequency in Hz
g           = 45.0; % antenna gain in dB
sigma       = 0.1; % radar cross section in m squared
b           = 5.0e+6; % radar operating bandwidth in Hz
nf          = 3.0; % noise figure in dB
loss        = 6.0; % radar losses in dB
range       = linspace(5e3,150e3,1000);%range : 5 km to 150 km
[pr, snr]   = radar_eq(pt, freq, g, sigma, b, nf, loss, range);
rangekm     = range ./ 1000;

subplot(3,1,3)
plot(rangekm,pr,'linewidth',1.5)
grid;
xlabel ('Detection range in km');
ylabel ('Pr in dBm');

%% 1. Active SAR (FM Transmitter) 
%FP01_CreateSatGeometry
%FP02_FindSwath

%% 2. Passive SAR (FM Receiver)
%FP01_CreateSatGeometry
%FP02_FindSwath

%% 3. Find GRP
%geodetic2aer

%% 4. Generate spatial sampling points (Targets)
%FP03_GenerateTargets

%% 5. Get ground reflectivity
%FP04_GetGroundReflect

%% 6. Define the source waveform
% this will be FM 

%% 7. Generate base chirp 
%sb = 

%% 8. Select the Testing value for testing the script
%Testing

%% 9. Approx azimuth of the satellite to calculate the antenna pattern
%azimuth

%% 10. Reference sqd that will be used for template match filtering
%FP05_CalcReflectionBi --> at GRP

%% 11. Generate the reflected signal from the entire swath
%FP05_CalcReflectionBi --> at TargetLat and TargetLon

%% 12. Plot the raw unfocused SAR signal



%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pr, snr] = radar_eq(pt, freq, g, sigma, b, nf, loss, range)
% 
% (SNR)o = (Pt*(G^2)*(λ^2)*σ)/((4π)^3)*k*T0*B*F*L*(R^4)
c               = 3.0e+8; % speed of light
lambda          = c / freq; % wavelength
p_peak          = 10*log10(pt); % convert peak power to dB
lambda_sqdb     = 10*log10(lambda^2); % compute wavelength square in dB
sigmadb         = 10*log10(sigma); % convert sigma to dB
four_pi_cub     = 10*log10((4.0 * pi)^3); % (4pi)^3 in dB
k_db            = 10*log10(1.38e-23); % Boltzmann's constant in dB
to_db           = 10*log10(290); % noise temp. in dB
b_db            = 10*log10(b); % bandwidth in dB
range_pwr4_db   = 10*log10(range.^4); % vector of target range^4 in dB
% Implement Equation (1.83)
num = p_peak + g + g + lambda_sqdb + sigmadb;% nominator
den = four_pi_cub + k_db + to_db + b_db + nf + loss + range_pwr4_db;% denominator
snr = num - den;
pr = num - (four_pi_cub + range_pwr4_db);
end