function [R,Ro,Idx,GRP] = F03_FindGRP(latSawthMid,lonSwathMid,Satlla,E)

% Create Range vector at the mid of the swath
R=[];

for etaIdx=1:length(latSawthMid)
    [~,~,R(etaIdx,:)] = geodetic2aer(latSawthMid(etaIdx),lonSwathMid(etaIdx),0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);
    %plot(R(etaIdx,:))
    %hold on
end

% Find eta with most symetrical pass
[~,etaIdx]=min(abs(R(:,end)-R(:,1)));
% Find GRP Range across the whole dwell time
R = R(etaIdx,:);
% Find Reference range at the centre of the dwell at ground refernece point(GRP) - Find Index of mid point of the dwell
[Ro,Idx] = min(R);
% Swath center point (Ground Reference Point GRP)
GRP = [latSawthMid(Idx),lonSwathMid(Idx),0]; 
end