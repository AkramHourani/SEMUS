R_SAR=Ro
varphi = acos( (-R_SAR.^2+(h+Re)^2+Re^2) ./ (2*(h+Re)*Re))
R_grd = varphi*Re
varphi = R_grd/Re

R_SAR = sqrt((h+Re)^2 + Re^2 - 2*(h+Re)*Re*cos(varphi))