function [R_SAR] = F07_SARRange(h,Re,R_grd)
varphi = R_grd/Re;
R_SAR = sqrt((h+Re)^2 + Re^2 - 2*(h+Re)*Re*cos(varphi));
end