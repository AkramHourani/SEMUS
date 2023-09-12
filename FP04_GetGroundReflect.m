function [a] = FP04_GetGroundReflect(Targetlat,Targetlon,latSawthL1,lonSwathL1,latSawthL2,lonSwathL2)
% To get the ground reflectivity based on the optical data 

% Extract the edges of the swath region
latEdge    = [min([latSawthL1;latSawthL2]), max([latSawthL1;latSawthL2])];
lonEdge    = [min([lonSwathL1;lonSwathL2]), max([lonSwathL1;lonSwathL2])];
FigOptical = figure('Position',[0 0 1000 1000]);                              % 1000 pixels
clf
%geoplot(Targetlat(:),Targetlon(:),'.'); 

%% Swath edge line

% Function : geolimits(latlim,lonlim) adjusts the limits of the current geographic axes, 
% geographic bubble chart, or map axes to include the latitude and longitude 
% limits specified by latlim and lonlim, respectively. 
% If the current axes is not a geographic axes, a geographic bubble chart, or a map axes, 
% or if there is no current axes, then the function creates a new geographic axes with the specified limits.
geolimits(latEdge,lonEdge)
geobasemap satellite
[latEdge, lonEdge] = geolimits; % Extract ground truth limits becuase the requested limits are not the same as the actual limits
ax = gca;
ax.Position=[0 0 1 1];
%pause(5);
disp("If the image is loaded press any key to continue...")
pause
Optical = rgb2gray(frame2im(getframe(ax)));
close (FigOptical)

%% Interpolating the optical map
latVec = linspace(latEdge(1),latEdge(2),size(Optical,1));
lonVec = linspace(lonEdge(1),lonEdge(2),size(Optical,2));
[lonVec,latVec] = meshgrid(lonVec,latVec);

% Function : B = flipud(A) returns A with its rows flipped in the up-down 
% direction (that is, about a horizontal axis).
latVec = flipud(latVec);

% Function : Vq = interp2(X,Y,V,Xq,Yq) returns interpolated values of 
% a function of two variables at specific query points using linear interpolation. 
a = double(interp2(lonVec,latVec,Optical,Targetlon,Targetlat,'nearest'));

end
