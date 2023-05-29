function MagArr = FP06_GetPSLR(a_az,Img,a_rg_0,a_rg_1,flag)

%% a. The border between mainlobes
idx_mg = a_rg_0:a_rg_1;

mag = 10*log10(Img(a_az,idx_mg));%%range
if flag ==1 
    mag = 10*log10(Img(idx_mg,a_az));%%azimuth
end
MagArr = zeros(1,18);

%% b. Finding the mainlobe
[PkAmp, PkTime] = findpeaks(mag,idx_mg);
[~,idx] = sort(PkAmp,'descend');
MagArr(1,2) = PkTime(idx(1));%mainlobe  
MagArr(1,5) = PkAmp(idx(1));%mainlobe 

%interpolation for plotting
int_id = idx_mg(1):0.001:idx_mg(end);
xSpline = interp1(idx_mg,mag,int_id,'spline');
[PkAmp, PkTime] = findpeaks(xSpline,int_id);
[~,idxp] = sort(PkAmp,'descend');
MagArr(1,8) = PkTime(idxp(1));%mainlobe  
MagArr(1,11) = PkAmp(idxp(1));%mainlobe

%% then we slice the signal with mainlobe as mid-value 
% c. Finding the first null (border for exclusion zone)

idx_mg_lf = a_rg_0:MagArr(1,2);

mag_lf = 10*log10(Img(a_az,idx_mg_lf));%%range
if flag ==1 
    mag_lf = 10*log10(Img(idx_mg_lf,a_az));%%azimuth
end

[NuAmp, NuTime] = findpeaks(max(mag_lf)-mag_lf,idx_mg_lf);

[~,idx] = sort(NuAmp,'ascend');
MagArr(1,1) = NuTime(idx(1));%first null 
MagArr(1,4) = max(mag)-NuAmp(idx(1));%left-side null 

%interpolation for plotting left-side null
int_id = idx_mg_lf(1):0.001:idx_mg_lf(end);
xSpline1 = interp1(idx_mg_lf,mag_lf,int_id,'spline');
[NuAmp, NuTime] = findpeaks(max(xSpline1)-xSpline1,int_id);
[~,idxp] = sort(NuAmp,'ascend');
MagArr(1,7) = NuTime(idxp(1));%first null  
MagArr(1,10) = max(xSpline1)-NuAmp(idxp(1));%left-side null 
 
%% d. Finding the second null (border for exclusion zone)
idx_mg_rg = MagArr(1,2):a_rg_1;
mag_rg = 10*log10(Img(a_az,idx_mg_rg));%%range
if flag ==1 
    mag_rg = 10*log10(Img(idx_mg_rg,a_az));%%azimuth
end

[NuAmp, NuTime] = findpeaks(max(mag_rg)-mag_rg,idx_mg_rg);
[~,idx] = sort(NuAmp,'ascend');
MagArr(1,3) = NuTime(idx(1));%mainlobe 
MagArr(1,6) = max(mag)-NuAmp(idx(1));%left-side null 

%interpolation for plotting right-side null
int_id = idx_mg_rg(1):0.001:idx_mg_rg(end);
xSpline1 = interp1(idx_mg_rg,mag_rg,int_id,'spline');
[NuAmp, NuTime] = findpeaks(max(xSpline1)-xSpline1,int_id);
[~,idxp] = sort(NuAmp,'ascend');
MagArr(1,9) = NuTime(idxp(1));%mainlobe  
MagArr(1,12) = max(xSpline1)-NuAmp(idxp(1));%left-side null 

%% to find nan mean
%MagArr(1,14) = nanmean(Img(a_az,MagArr(1,3):a_rg_1));%%range

img_comb = [Img(a_az,a_rg_0:MagArr(1,1)) Img(a_az,MagArr(1,3):a_rg_1)];

if flag ==1 
    img_comb = [Img(a_rg_0:MagArr(1,1),a_az);Img(MagArr(1,3):a_rg_1,a_az)];
end
MagArr(1,13) = nanmean(img_comb);%%range
%p_sl = 10*log10((MagArr(1,13)+MagArr(1,14))/2);
p_sl = 10*log10(MagArr(1,13));
MagArr(1,15) = p_sl;


pslr = MagArr(1,5) - MagArr(1,15);
MagArr(1,16) = pslr;
MagArr(1,17) = 10*log10(max(img_comb));
MagArr(1,18) = 10*log10(min(img_comb));







end