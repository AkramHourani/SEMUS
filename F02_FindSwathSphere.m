function [latSawthMid,lonSwathMid,Swathwidth,latSawthL1,lonSwathL1,latSawthL2,lonSwathL2,R_grd_o]=F02_FindSwathV2(Satlla,RadPar,Re,h)

Incident_o = asind((Re+h)/Re*sind(RadPar.AntOffNadir));
varphi_o = Incident_o-RadPar.AntOffNadir;
R_grd_o = varphi_o*Re/180*pi;

Incident_edge = asind((Re+h)/Re*sind(RadPar.AntOffNadir+RadPar.SwathWidthDeg /2));
varphi_edge = Incident_edge-(RadPar.AntOffNadir+RadPar.SwathWidthDeg /2);
R_grd_edge = varphi_edge*Re/180*pi;

Swathwidth = abs((R_grd_edge-R_grd_o)*2);

% This is to compute the approximate azimuth of the swath (also the
% satellite azimuth -> i.e. the direction of motion of the satellite
if RadPar.Left == 0
    az=90;
else
    az=-90;
end

[latSawthMid,lonSwathMid] = reckon(Satlla(:,1),Satlla(:,2),varphi_o,az);
[latSawthL1,lonSwathL1]  = reckon(latSawthMid,lonSwathMid,Swathwidth/2/Re*180/pi,90);
[latSawthL2,lonSwathL2]  = reckon(latSawthMid,lonSwathMid,Swathwidth/2/Re*180/pi,-90);

end