function [sigma] = F05_GetGroundReflect(Targetlat,Targetlon,latSwathL1,lonSwathL1,latSwathL2,lonSwathL2)
% Here we get the gournd reflectivity based on the optical data 
% Extract the edges of the swath region
latEdge = [min([latSwathL1;latSwathL2]), max([latSwathL1;latSwathL2])];
lonEdge = [min([lonSwathL1;lonSwathL2]), max([lonSwathL1;lonSwathL2])];
FigOptical=figure('Position',[0 0 1000 1000]);         % 1000 pixels
clf
%geoplot(Targetlat(:),Targetlon(:),'.');               % Swath edge line
geolimits(latEdge,lonEdge)
geobasemap satellite
% Extract ground truth limits becuase the requested limits are not the same as the actual limits
[latEdge, lonEdge] = geolimits;                         
ax=gca;
ax.Position=[0 0 1 1];
%pause(5);
disp("If the image is loaded press any key to continue...")
pause
Optical=rgb2gray(frame2im(getframe(ax)));
close (FigOptical)
%% Interpolating the optical map
latVec = linspace(latEdge(1),latEdge(2),size(Optical,1));
lonVec = linspace(lonEdge(1),lonEdge(2),size(Optical,2));
[lonVec,latVec]=meshgrid(lonVec,latVec);
latVec = flipud(latVec);
sigma = double(interp2(lonVec,latVec,Optical,Targetlon,Targetlat,'nearest'));