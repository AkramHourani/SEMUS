clc; clear; close all hidden
load('Points.mat')
%load('Sydney.mat')
%load('SAR_Image1P.mat')

%% This is a row-wise FFT / IFFT
fft1d2  = @ (x) fftshift(fft(fftshift(x,2),[],2),2);
ifft1d2 = @ (x) ifftshift(ifft(ifftshift(x,2),[],2),2);
fft1d1  = @ (x) fftshift(fft(fftshift(x,1),[],1),1);
ifft1d1 = @ (x) ifftshift(ifft(ifftshift(x,1),[],1),1);

%% Plotting raw time domain signal
figure_P1a

%% Step 1: Range Compression
So = fft1d2(sqd); % FFT the time domain signal (FFT along each eta row)

% This is the template for the Match Filter
tau     = 0;
sb      = exp(1j*pi *   (-2*RadPar.fo * tau + RadPar.K*(FastTime-tau).^2   )    ) ...
            .*(FastTime>(-RadPar.T/2+tau)).*(FastTime<(RadPar.T/2+tau));
G       = conj(fftshift(fft(fftshift(sb))));
Src     = repmat(G,size(So,1),1).*So; %Equation 5.5
src     = ifft1d2(Src); %Equation 5.6

%% Plotting the Real component of the range matched filter
nexttile;figure_P1b

%% Generate the reference signal based on the ground reference point
So      = fft1d2(sqd_ref); % FFT the time domain signal (FFT along each eta row)
% This is the template for the Match Filter : complex conjugate of the
% input signal
tau     = 0;
sb      = exp(1j*pi *   (-2*RadPar.fo * tau + RadPar.K*(FastTime-tau).^2   )    ) ...
    .*(FastTime>(-RadPar.T/2+tau)).*(FastTime<(RadPar.T/2+tau));
G       = conj(fftshift(fft(fftshift(sb))));  % G(f)
Src     = repmat(G,size(So,1),1).*So; %Equation 5.5
src_ref = ifft1d2(Src); %Equation 5.6 

%% Plotting the range-compressed image
nexttile;figure_P1c

%% Step 2 Azimuth FFT
S2_ref  = fft1d1(src_ref);
S2      = fft1d1(src);

nexttile;figure_P1d

%% Step 3 Range cell migration compensation
DeltaR = RSoI+RI-Ro; % The migration of the GRP

nexttile;figure_P1e

% shifting range cells (Needs to be update for Bi Static SAR)
RangeBin    = RadPar.ts*c/2;
NbinsShift  = -round(DeltaR/RangeBin)*2;
% RangeBin  = RadPar.ts*c;
% NbinsShift = -round(DeltaR/RangeBin);
for eta=1:etaTotal
    S2(eta,:)       = circshift(S2(eta,:),NbinsShift(eta));
    S2_ref(eta,:)   = circshift(S2_ref(eta,:),NbinsShift(eta));
end

nexttile;figure_P1f

%% Step 4/5 Azimuth compression
%%% Empirical Method %%%%
midpoint    = round(size(S2_ref,2)/2)+1;
S2_ref      = repmat(S2_ref(:,midpoint),1,size(S2,2));
S3          = S2 .* conj(S2_ref);

%% Step 5 Azimuth IFFT
sSLC        = ifft1d1(S3);

%% Plotting Compressed Image
%CalibrationR = 2;
%CalibrationAz = 160;
%Img = abs(sSLC)./max(abs(sSLC),[],"all");
Range       = (-(numel(FastTime)/2)*RangeBin:RangeBin:(numel(FastTime)/2-1)*RangeBin);
[~,El,~]    = geodetic2aer(GRP(1),GRP(2),0,SatllaSoI(1,1),SatllaSoI(1,2),SatllaSoI(1,3),E); %% I changed Satlla to SatllaSoI
RangeGround = Range/cosd(El);
speed       = mean(sqrt(sum((diff(SatECISoI,[],2)).^2)) /Param.ts);%% I changed SatECI to SatECISoI
CrossRange  = (1:etaTotal)*Param.ts*speed/1000;

J = abs(sSLC);
J = J./max(J,[],"all");
J = imresize(J,2);
J = imgaussfilt(J ,2); % Smoothing filter
J = imadjust(J,[0 0.5]);

figure_9b
% figure_P7
% create_vid_4

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

% %%
% % Transform to local cartesian coordinates