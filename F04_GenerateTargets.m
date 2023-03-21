function [Targetlat,Targetlon]= F04_GenerateTargets(latSawthL1,lonSwathL1,latSawthL2,lonSwathL2,Param)
% This loop will walk step by step along the track to generate the sampling
% points (arc) between the two edges based on the great cricle

Targetlat=zeros(length(latSawthL1),Param.NtargetsRange+1);
Targetlon=zeros(length(lonSwathL1),Param.NtargetsRange+1);
etaTotal=length(latSawthL1);                    % Eta - Azimuth 

for eta=1:etaTotal

% Option 1: random sampling
% Targetlat(eta,:)=rand(1,Param.NtargetsRange+1)*(latSawthL2(eta)-latSawthL1(eta))+latSawthL1(eta);
% Targetlon(eta,:)=rand(1,Param.NtargetsRange+1)*(lonSwathL2(eta)-lonSwathL1(eta))+lonSwathL1(eta);

% Option 2: regular sampling creats some visulaizaiton issues
[Targetlat(eta,:),Targetlon(eta,:)] = gcwaypts(latSawthL1(eta),lonSwathL1(eta),latSawthL2(eta),lonSwathL2(eta),Param.NtargetsRange); 

% Option 3: perturbed grid samples
% [Targetlat(eta,:),Targetlon(eta,:)] = gcwaypts(latSawthL1(eta),lonSwathL1(eta),latSawthL2(eta),lonSwathL2(eta),Param.NtargetsRange); 
% LatStep = abs(mean(diff(Targetlat(eta,:)))); % This is the spacing between samples
% LonStep = abs(mean(diff(Targetlon(eta,:))));
% 
% % Here we add the perturbation
% Targetlat(eta,:) = rand(size(Targetlat(eta,:))).*LatStep+Targetlat(eta,:);
% Targetlon(eta,:) = rand(size(Targetlon(eta,:))).*LonStep+Targetlon(eta,:);
end
end