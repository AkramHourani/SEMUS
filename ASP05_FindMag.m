%%startfromhere%%%%%%%%%%%%%%%%%%%%
close all;
% azres = 3.0231; %calculated from AS04
% slantrngres = 9.9931;%calculated from AS04
%% 1. Test update Using fast peak find

figure
%1.a Load the input image and determine the size
Imgslc = abs(sSLC);
[br,bc] = size(Imgslc);

%1.b Determine the threshold value
pix_value = 1e-14;

%1.c Add pixel for surrounding image
Img = zeros(br+20,bc+20);
Img(:,1:10) = pix_value;%add pixel left	
Img(1:10,:) = pix_value;%add pixel top
Img(:,bc+11:bc+20) = pix_value;%add pixel right
Img(br+11:br+20,:) = pix_value;%add pixel bottom	
Img(11:br+10,11:bc+10) = Imgslc(:,:);%add pixel bottom	

%1.d Collect all detected peaks for all region
% the peaks are collected in 'p' where [x,y] = [p(odd),p(even)]
p = FP08_FastPeakFind(Img);
imagesc(Img); 
hold on
plot(p(1:2:end),p(2:2:end),'y+')
hold on
% save the x and y position of deected peaks
x_axs = p(1:2:end);%range, odd numbered is x axis
y_axs = p(2:2:end);%azimuth, even numbered is y axis

% create an array : x,y with z(amplitude)
val_main    = zeros(length(p)/2,3);
for ty = 1 : length(p)/2
	val_main(ty,1) = x_axs(ty);
	val_main(ty,2) = y_axs(ty);
	val_main(ty,3) = Img(y_axs(ty), x_axs(ty));
end
% sort the array, based on the z(amplitude) from the max to min 
val_main = sortrows(val_main,3,'descend');

%% 2. Create smaller border to find local peak (the intersection between horizontal and vertical line of peaks) 
%2.a. Find vertical border based on X axis
% histogram is used to find peak occurences
figure
hx = histogram(x_axs,round(bc/Param.NtargetsRange));
title('hx');
hbinx(:,1) = round(hx.BinEdges);
hbinx(:,2) = [0 hx.BinCounts];

brds = 2.5;   %border distance
lim_az      = etaTotal;
lim_rg      = length(FastTime);

hbin2x = sortrows(hbinx,2,'descend');% sorting
hbin3x = hbin2x(1:Param.NtargetsRange,:);% get the first 11 large values
hbin4x = sortrows(hbin3x,1,'ascend');% sort again based on bins

%2.b Create cell border
%2.b.1 For the first row
hbin4x(1,3) = hbin4x(1,1) - round(minus(hbin4x(1,1),0)/brds);
hbin4x(1,4) = hbin4x(1,1) + round(minus(hbin4x(2,1),hbin4x(1,1))/brds);

%2.b.2 for the last row
hbin4x(length(hbin4x),3) = hbin4x(length(hbin4x),1) -  round(minus(hbin4x(length(hbin4x),1),hbin4x(length(hbin4x)-1,1))/brds);
hbin4x(length(hbin4x),4) = hbin4x(length(hbin4x),1) +  round(minus(lim_rg,hbin4x(length(hbin4x),1))/brds);

%2.b.3 for rows other than first and last
for tg = 2 : length(hbin4x)-1 
	hbin4x(tg,3) = hbin4x(tg,1) - round(minus(hbin4x(tg,1),hbin4x(tg-1,1))/brds);
	hbin4x(tg,4) = hbin4x(tg,1) + round(minus(hbin4x(tg+1,1),hbin4x(tg,1))/brds);
end 

%% 2.c. Find horizontal border based on Y axis
% histogram is used to find peak occurences
figure
hy = histogram(y_axs,round(2*br/Param.NtargetsAz));
title('hy');
hbiny(:,1) = round(hy.BinEdges);
hbiny(1,1) = 1;%can not zero
hbiny(:,2) = [hy.BinCounts 0];%count of peaks detected

for ht = 1:length(hbiny)    
    if ht == length(hbiny)%reaches the end of rows
        hbiny(ht,4) = 1;%equals to 1
    else
        hbiny(ht,4) = hbiny(ht,2)/hbiny(ht+1,2);
    end
    
    if ht == 1%the beginning, the first row
        hbiny(ht,3) = hbiny(ht,2)/hbiny(ht,2);%equals to 1
    else
        hbiny(ht,3) = hbiny(ht,2)/hbiny(ht-1,2);
    end
end 

hbiny(isinf(hbiny)|isnan(hbiny)) = 0;
hbin2y = sortrows(hbiny,2,'descend');% sorting

% get first 11 large values
hbin3y = hbin2y(1:round(1.5*Param.NtargetsAz),:);%extend the search to 1.5 time targets
hbin3y = sortrows(hbin3y,1,'ascend');
for ht = 1:length(hbin3y)    
    
    if ht == length(hbin3y)%reaches the end of rows
        hbin3y(ht,5) = hbin3y(ht,1);
    else
        hbin3y(ht,5) = abs(hbin3y(ht,1)-hbin3y(ht+1,1));
    end
end 

hbin3y = sortrows(hbin3y,5,'descend');

%sort again based on bins
hbin4y = hbin3y(hbin3y(:,5)> round(bc/Param.NtargetsRange));
hbin4y = sortrows(hbin4y,1,'ascend');

%2.d Create cell border
%2.d.1 For the first row
hbin4y(1,3) = hbin4y(1,1) - round(minus(hbin4y(1,1),0)/brds);%limit before
hbin4y(1,4) = hbin4y(1,1) + round(minus(hbin4y(2,1),hbin4y(1,1))/brds);%limit after

%2.d.2 for the last row
hbin4y(length(hbin4x),3) =  hbin4y(length(hbin4y),1) -  round(minus(hbin4y(length(hbin4y),1),hbin4y(length(hbin4y)-1,1))/brds);%limit before
hbin4y(length(hbin4y),4) = br+20;

%2.d.3 for rows other than first and last
for tg = 2 : length(hbin4y)-1 
	hbin4y(tg,3) = hbin4y(tg,1) - round(minus(hbin4y(tg,1),hbin4y(tg-1,1))/brds);
	hbin4y(tg,4) = hbin4y(tg,1) + round(minus(hbin4y(tg+1,1),hbin4y(tg,1))/brds);
end 

%% Plotting cells
Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);

imagesc(Img); 
%colormap sky
%title("Target Points",'interpreter','latex');
xlabel('Fast time','interpreter','Tex')
ylabel('Slow time','interpreter','Tex');

nd = 1;
ne = 1;
trg_val = zeros;

for gy = 1:11 
	for gx = 1:11 
		hold on
       
        Imgs = Img(hbin4y(gy,3):hbin4y(gy,4),hbin4x(gx,3):hbin4x(gx,4)); 
	    ptest = FP08_FastPeakFind(Imgs);
	    x_cell = hbin4x(gx,3) + ptest(1:2:end);%range, odd numbered is x axis
	    y_cell = hbin4y(gy,3) + ptest(2:2:end);%azimuth, even numbered is y axis

        val_cell = zeros(length(ptest)/2,3);
		
        for ty = 1 : length(ptest)/2
			val_cell(ty,1:3) = [x_cell(ty) y_cell(ty) Img(y_cell(ty), x_cell(ty))];
		end

	    %now find max
	    [mx,ix] = max(val_cell(:,3));
        trg_val(nd,1:5) = [gx gy val_cell(ix,:)];%save results to array

        %%plotting 
        %plot small cell borders
 		rt = rectangle('Position',[hbin4x(gx,3) hbin4y(gy,3) minus(hbin4x(gx,4),hbin4x(gx,3)) minus(hbin4y(gy,4),hbin4y(gy,3))]);
        rt.EdgeColor = 'y';
        
	    %plot peaks inside the small cells
        plot(x_cell  ,y_cell,'y+');
        
        %plot plus sign at target points
        str = sprintf('[%d,%d]', gx, gy);
	    plot(val_cell(ix,1),val_cell(ix,2) ,'r+');
        
        %put label text
        text(val_cell(ix,1)-10, val_cell(ix,2)-20, str, 'FontSize', 7, 'Color', 'w');
        %%%%%%%%%%%

        nd = nd + 1;
    end
    nd = nd+1;
end

%% Target Point at GRP and Calculating PSLR

%1. Point A : start scan, range compressed
a_pt = [6 2];%range, azimuth
[a_rg,a_az,a_rg_0,a_rg_1,a_az_0,a_az_1] = FP09_PlotMag(a_pt(1),a_pt(2),trg_val,length(FastTime),etaTotal);
MagArr(1,:) = FP06_GetPSLR(a_az,Img,a_rg_0,a_rg_1,0);

%2. GRP : find the grp, range compressed, azimuth compressed
grp_pt = [6 6];
[grp_rg,grp_az,grp_rg_0,grp_rg_1,grp_az_0,grp_az_1] = FP09_PlotMag(grp_pt(1),grp_pt(2),trg_val,length(FastTime),etaTotal);
MagArr(2,:) = FP06_GetPSLR(grp_az,Img,grp_rg_0,grp_rg_1,0);
MagArr(5,:) = FP06_GetPSLR(grp_rg,Img,grp_az_0,grp_az_1,1);

%3. Point B : finish scan, range compressed
b_pt = [6 8];
[b_rg,b_az,b_rg_0,b_rg_1,b_az_0,b_az_1] = FP09_PlotMag(b_pt(1),b_pt(2),trg_val,length(FastTime),etaTotal);
MagArr(3,:) = FP06_GetPSLR(b_az,Img,b_rg_0,b_rg_1,0);

%4. Point C : near edge, azimuth compressed
c_pt = [2 6];
[c_rg,c_az,c_rg_0,c_rg_1,c_az_0,c_az_1] = FP09_PlotMag(c_pt(1),c_pt(2),trg_val,length(FastTime),etaTotal);
MagArr(4,:) = FP06_GetPSLR(c_rg,Img,c_az_0,c_az_1,1);

%5. Point D : finish scan, azimuth compressed
d_pt = [10 6];
[d_rg,d_az,d_rg_0,d_rg_1,d_az_0,d_az_1] = FP09_PlotMag(d_pt(1),d_pt(2),trg_val,length(FastTime),etaTotal);
MagArr(6,:) = FP06_GetPSLR(d_rg,Img,d_az_0,d_az_1,1);

set(gca,'LooseInset',get(gca,'TightInset'));
print(h_Fig, '-dpng','-r600','Figure13')
movefile('Figure13.png','Figures')

%% Create Plot : Range Compression
Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');

nexttile %subplot(3,2,1)
FP07_LinePlot(1,1,Img,a_az,1:length(FastTime),MagArr(1,8),MagArr(1,11),sprintf('$L_1$'),0,0,0,0,0,0,0,[0.5, 0.1, 0])

nexttile %subplot(3,2,2)
FP07_LinePlot(2,1,Img,a_az,a_rg_0:a_rg_1,MagArr(1,8),MagArr(1,11),sprintf('$P_1$ PSLR = %.1f dB',MagArr(1,16)),0,0, MagArr(1,7), MagArr(1,9), MagArr(1,15), MagArr(1,17), MagArr(1,18),[0.5, 0.1, 0])

nexttile %subplot(3,2,3)
FP07_LinePlot(1,1,Img,grp_az,1:length(FastTime),MagArr(2,8),MagArr(2,11),sprintf('$LGRP_{r}$'),0,1,0,0,0,0,0,[0.45, 0.1, 0])

nexttile %subplot(3,2,4)
FP07_LinePlot(2,1,Img,grp_az,grp_rg_0:grp_rg_1,MagArr(2,8),MagArr(2,11),sprintf('$GRP_{r}$ PSLR = %.1f dB',MagArr(2,16)),0,0, MagArr(2,7), MagArr(2,9), MagArr(2,15), MagArr(2,17), MagArr(2,18),[0.5, 0.1, 0])

nexttile %subplot(3,2,5)
FP07_LinePlot(1,1,Img,b_az,1:length(FastTime),MagArr(3,8),MagArr(3,11),sprintf('$L_2$'),1,0,0,0,0,0,0,[0.5, 0.1, 0])

nexttile %subplot(3,2,6)
FP07_LinePlot(2,1,Img,b_az,b_rg_0:b_rg_1,MagArr(3,8),MagArr(3,11),sprintf('$P_2$ PSLR = %.1f dB',MagArr(3,16)),1,0, MagArr(3,7), MagArr(3,9), MagArr(3,15), MagArr(3,17), MagArr(3,18),[0.5, 0.1, 0])

print(h_Fig, '-dpng','-r600','Figure11')
movefile('Figure11.png','Figures')

%% Create Plot Azimuth Compression
Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');

nexttile %subplot(3,2,1)
FP071_LinePlot(1,    2,    Img, c_rg, 1:etaTotal,MagArr(4,8),MagArr(4,11),sprintf('$L_3$'),0,0,0,0,0,0,0,[0.5, 0.1, 0])

nexttile %subplot(3,2,2)
FP071_LinePlot(2,2,Img,c_rg,c_az_0:c_az_1,MagArr(4,8),MagArr(4,11),sprintf('$P_3$'),0,0, MagArr(4,7), MagArr(4,9), MagArr(4,15), MagArr(4,17), MagArr(4,18),[0.5, 0.1, 0])

nexttile %subplot(3,2,3)
FP07_LinePlot(1,2,Img,grp_rg,1:etaTotal,MagArr(5,8),MagArr(5,11),sprintf('$LGRP_a$'),0,1,0,0,0,0,0,[0.45, 0.1, 0])

nexttile %subplot(3,2,4)
FP071_LinePlot(2,2,Img,grp_rg,grp_az_0:grp_az_1,MagArr(5,8),MagArr(5,11),sprintf('$GRP_{a}$'),0,0, MagArr(5,7), MagArr(5,9), MagArr(5,15), MagArr(5,17), MagArr(5,18),[0.45, 0.1, 0])

nexttile %subplot(3,2,5)
FP071_LinePlot(1,2,Img,d_rg,1:etaTotal,MagArr(6,8),MagArr(6,11),sprintf('$L_4$'),2,0,0,0,0,0,0,[0.5, 0.1, 0])

nexttile %subplot(3,2,6)
FP071_LinePlot(2,2,Img,d_rg,d_az_0:d_az_1,MagArr(6,8),MagArr(6,11),sprintf('$P_4$'),2,0, MagArr(6,7), MagArr(6,9), MagArr(6,15), MagArr(6,17), MagArr(6,18),[0.5, 0.1, 0])

print(h_Fig, '-dpng','-r600','Figure10')
movefile('Figure10.png','Figures')
