%%startfromhere%%%%%%%%%%%%%%%%%%%%
close all;
% azres = 3.0231; %calculated from AS04
% slantrngres = 9.9931;%calculated from AS04
%% 1. using fast peak find
figure
Imgslc = abs(sSLC);

[br,bc] = size(Imgslc);
pix_value = mean(Imgslc(:)); %1e-14;

%add pixel for surrounding image
Img = zeros(br+20,bc+20);
Img(:,1:10) = pix_value;%add pixel left	
Img(1:10,:) = pix_value;%add pixel top
Img(:,bc+11:bc+20) = pix_value;%add pixel right
Img(br+11:br+20,:) = pix_value;%add pixel bottom	
Img(11:br+10,11:bc+10) = Imgslc(:,:);%add pixel bottom	


lim_az = etaTotal;
lim_rg = length(FastTime);

p = FP08_FastPeakFind(Img);
imagesc(Img); 
hold on

x_axs = p(1:2:end);%range, odd numbered is x axis
y_axs = p(2:2:end);%azimuth, even numbered is y axis

val_main = zeros(length(p)/2,3);

for ty = 1 : length(p)/2
	val_main(ty,1) = x_axs(ty);
	val_main(ty,2) = y_axs(ty);
	val_main(ty,3) = Img(y_axs(ty), x_axs(ty));
end

val_main = sortrows(val_main,3,'descend');

brds = 2.5;   %border distance

plot(p(1:2:end),p(2:2:end),'y+')
hold on
%% X axis
figure
hx = histogram(x_axs,round(bc/Param.NtargetsRange));
title('hx');
hbinx(:,1) = round(hx.BinEdges);
hbinx(:,2) = [0 hx.BinCounts];

% sorting
hbin2x = sortrows(hbinx,2,'descend');

% get first 11 large values
hbin3x = hbin2x(1:11,:);

%sort again based on bins
hbin4x = sortrows(hbin3x,1,'ascend');

%create cell border
%for the first row
hbin4x(1,3) = hbin4x(1,1) - round(minus(hbin4x(1,1),0)/brds);
hbin4x(1,4) = hbin4x(1,1) + round(minus(hbin4x(2,1),hbin4x(1,1))/brds);

%for the last row
hbin4x(length(hbin4x),3) = hbin4x(length(hbin4x),1) -  round(minus(hbin4x(length(hbin4x),1),hbin4x(length(hbin4x)-1,1))/brds);
hbin4x(length(hbin4x),4) = hbin4x(length(hbin4x),1) +  round(minus(lim_rg,hbin4x(length(hbin4x),1))/brds);
%hbin4x(length(hbin4x),4) = length(FastTime);

%for rows other than first and last
for tg = 2 : length(hbin4x)-1 
	hbin4x(tg,3) = hbin4x(tg,1) - round(minus(hbin4x(tg,1),hbin4x(tg-1,1))/brds);
	hbin4x(tg,4) = hbin4x(tg,1) + round(minus(hbin4x(tg+1,1),hbin4x(tg,1))/brds);
end 

%% Y axis
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
    
    if ht == 1
        hbiny(ht,3) = hbiny(ht,2)/hbiny(ht,2);%equals to 1
    else
        hbiny(ht,3) = hbiny(ht,2)/hbiny(ht-1,2);
    end
end 

hbiny(isinf(hbiny)|isnan(hbiny)) = 0;

% sorting
hbin2y = sortrows(hbiny,2,'descend');

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
%hbin4y = hbin3y(:,5) > round(bc/Param.NtargetsRange) ;
hbin4y = hbin3y(hbin3y(:,5)> round(bc/Param.NtargetsRange));
hbin4y = sortrows(hbin4y,1,'ascend');

%create cell border

%for the first row
hbin4y(1,3) = hbin4y(1,1) - round(minus(hbin4y(1,1),0)/brds);%limit before
hbin4y(1,4) = hbin4y(1,1) + round(minus(hbin4y(2,1),hbin4y(1,1))/brds);%limit after

%for the last row
hbin4y(length(hbin4x),3) =  hbin4y(length(hbin4y),1) -  round(minus(hbin4y(length(hbin4y),1),hbin4y(length(hbin4y)-1,1))/brds);%limit before
%hbin4y(length(hbin4y),4) =  hbin4y(length(hbin4y),1) +  round(minus(lim_az,hbin4y(length(hbin4y),1))/brds);%limit after
hbin4y(length(hbin4y),4) = br+20;

%for rows other than first and last
for tg = 2 : length(hbin4y)-1 
	hbin4y(tg,3) = hbin4y(tg,1) - round(minus(hbin4y(tg,1),hbin4y(tg-1,1))/brds);
	hbin4y(tg,4) = hbin4y(tg,1) + round(minus(hbin4y(tg+1,1),hbin4y(tg,1))/brds);
end 

%%

%%%%creating cells
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
			val_cell(ty,1) = x_cell(ty);
			val_cell(ty,2) = y_cell(ty);
			val_cell(ty,3) = Img(y_cell(ty), x_cell(ty));
		end

	    %now find max
	    [mx,ix] = max(val_cell(:,3));

        trg_val(nd,1) = gx;
        trg_val(nd,2) = gy;
        trg_val(nd,3) = val_cell(ix,1);
        trg_val(nd,4) = val_cell(ix,2);
        trg_val(nd,5) = val_cell(ix,3);
        
        %%plotting 
        %plot small cell borders
 		%rt = rectangle('Position',[hbin4x(gx,3) hbin4y(gy,3) minus(hbin4x(gx,4),hbin4x(gx,3)) minus(hbin4y(gy,4),hbin4y(gy,3))]);
        %rt.EdgeColor = 'y';
        
	    %plot peaks inside the small cells
        %plot(x_cell  ,y_cell,'y+');
        
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



%% Target Point at GRP

%%Range Compressed

%1. Point A : start scan, range compressed
a_rpt = 6;
a_zpt = 2;

[a_rg,a_az,a_rg_0,a_rg_1,a_az_0,a_az_1] = FP09_PlotMag(a_rpt,a_zpt,trg_val,length(FastTime),etaTotal);

MagArr(1,:) = FP06_GetPSLR(a_az,Img,a_rg_0,a_rg_1,0);

%% d. Calculating PSLR

%2. GRP : find the grp, range compressed, azimuth compressed
rg_pt = 6;
az_pt = 6;

[grp_rg,grp_az,grp_rg_0,grp_rg_1,grp_az_0,grp_az_1] = FP09_PlotMag(rg_pt,az_pt,trg_val,length(FastTime),etaTotal);

MagArr(2,:) = FP06_GetPSLR(grp_az,Img,grp_rg_0,grp_rg_1,0);
MagArr(5,:) = FP06_GetPSLR(grp_rg,Img,grp_az_0,grp_az_1,1);

%3. Point B : finish scan, range compressed
b_rpt = 6;
b_zpt = 8;

[b_rg,b_az,b_rg_0,b_rg_1,b_az_0,b_az_1] = FP09_PlotMag(b_rpt,b_zpt,trg_val,length(FastTime),etaTotal);

MagArr(3,:) = FP06_GetPSLR(b_az,Img,b_rg_0,b_rg_1,0);

%4. Point C : near edge, azimuth compressed
c_rpt = 2;
c_zpt = 6;

[c_rg,c_az,c_rg_0,c_rg_1,c_az_0,c_az_1] = FP09_PlotMag(c_rpt,c_zpt,trg_val,length(FastTime),etaTotal);

MagArr(4,:) = FP06_GetPSLR(c_rg,Img,c_az_0,c_az_1,1);

%5. Point D : finish scan, azimuth compressed
d_rpt = 10;
d_zpt = 6;

[d_rg,d_az,d_rg_0,d_rg_1,d_az_0,d_az_1] = FP09_PlotMag(d_rpt,d_zpt,trg_val,length(FastTime),etaTotal);

MagArr(6,:) = FP06_GetPSLR(d_rg,Img,d_az_0,d_az_1,1);

%start draw lines
% yline(a_az,'-.','LineWidth', 2,'Color',"#D95319");%L1
% yline(grp_az,'-.','LineWidth', 2,'Color',"b");%LGRPa
% yline(b_az,'-.','LineWidth', 2, 'Color',"#D95319");%L2
% xline(c_rg,'-.','LineWidth', 2,'Color',"#006400");%L3
% xline(grp_rg,'-.','LineWidth', 2,'Color',"b");%LGRPr
% xline(d_rg,'-.','LineWidth', 2,'Color',"#006400");%L4
% 
% plot(a_rg,a_az,'kO','MarkerSize',10)%P1
% plot(b_rg,b_az,'kO','MarkerSize',10)
% plot(grp_rg,grp_az,'kO','MarkerSize',10)%GRP
% plot(c_rg,c_az,'kO','MarkerSize',10)%P3
% plot(d_rg,d_az,'kO','MarkerSize',10)%P4
% 
% text(grp_rg-100, grp_az-30,'GRP', 'FontSize', 13, 'Color', 'k','Interpreter','latex');
% text(grp_rg+10, grp_az-280,'$LGRP_{a}$', 'FontSize', 13, 'Color', 'b','Interpreter','latex');
% text(grp_rg+180, grp_az-40,'$LGRP_{r}$', 'FontSize', 13, 'Color', 'b','Interpreter','latex');
% text(c_rg+10, a_az+40,'$L_{1}$', 'FontSize', 13, 'Color', "#D95319", 'Interpreter','latex');
% text(c_rg+10, b_az-40,'$L_{2}$', 'FontSize', 13, 'Color', "#D95319", 'Interpreter','latex');
% text(c_rg+8,  grp_az-180,'$L_{3}$', 'FontSize', 13, 'Color', '#006400','Interpreter','latex');
% text(d_rg-60,  grp_az-180,'$L_{4}$', 'FontSize', 13, 'Color', '#006400','Interpreter','latex');
% text(a_rg+20, a_az+40,'$P_{1}$', 'FontSize', 13, 'Color', 'k', 'Interpreter','latex');
% text(a_rg+20, b_az+40,'$P_{2}$',  'FontSize', 13, 'Color', 'k', 'Interpreter','latex');
% text(c_rg+20, c_az+30,'$P_{3}$', 'FontSize', 13, 'Color', 'k','Interpreter','latex');
% text(d_rg+20, c_az+30,'$P_{4}$', 'FontSize', 13, 'Color', 'k','Interpreter','latex');
%%draw lines ends here

set(gca,'LooseInset',get(gca,'TightInset'));
print(h_Fig, '-dpng','-r600','Figure13')
movefile('Figure13.png','Figures')

%% Find Magnitude 1-dimension
Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
t = tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');

nexttile %subplot(3,2,1)
idx_mg = 1:length(FastTime);
mag = 10*log10(Img(a_az,idx_mg));
subtt = sprintf('$L_1$');
FP07_LinePlot(1,idx_mg,mag,MagArr(1,8),MagArr(1,11),subtt,0,0)

title(subtt, 'Units', 'normalized', 'Position', [0.5, 0.1, 0],'horizontalAlignment', 'left','interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','k');

nexttile %subplot(3,2,2)
idx_mg = a_rg_0:a_rg_1;
mag = 10*log10(Img(a_az,idx_mg));

FP07_LinePlot(2,idx_mg,mag,MagArr(1,8),MagArr(1,11),subtt,0,0, MagArr(1,7), MagArr(1,9), MagArr(1,15), MagArr(1,17), MagArr(1,18))
subtt = sprintf('$P_1$ PSLR = %.1f dB',MagArr(1,16));
title(subtt, 'Units', 'normalized', 'Position', [0.1, 0.1, 0],'horizontalAlignment', 'left','interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','k');

nexttile %subplot(3,2,3)
idx_mg = 1:length(FastTime);
mag = 10*log10(Img(grp_az,idx_mg));

FP07_LinePlot(1,idx_mg,mag,MagArr(2,8),MagArr(2,11),subtt,0,1)
subtt = sprintf('$LGRP_{r}$');
title(subtt, 'Units', 'normalized', 'Position', [0.45, 0.1, 0],'horizontalAlignment', 'left','interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','k');

nexttile %subplot(3,2,4)
idx_mg = grp_rg_0:grp_rg_1;
mag = 10*log10(Img(grp_az,idx_mg));

%MagArr(2,7) =  MagArr(2,7)+2;
%MagArr(2,3) =  MagArr(2,2)+2;
FP07_LinePlot(2,idx_mg,mag,MagArr(2,8),MagArr(2,11),subtt,0,0, MagArr(2,7), MagArr(2,9), MagArr(2,15), MagArr(2,17), MagArr(2,18))
subtt = sprintf('$GRP_{r}$ PSLR = %.1f dB',MagArr(2,16));
title(subtt, 'Units', 'normalized', 'Position', [0.1, 0.1, 0],'horizontalAlignment', 'left','interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','k');

nexttile %subplot(3,2,5)
idx_mg = 1:length(FastTime);
mag = 10*log10(Img(b_az,idx_mg));

FP07_LinePlot(1,idx_mg,mag,MagArr(3,8),MagArr(3,11),subtt,1,0)
subtt = sprintf('$L_2$');
title(subtt, 'Units', 'normalized', 'Position', [0.5, 0.1, 0],'horizontalAlignment', 'left','interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','k');

nexttile %subplot(3,2,6)
idx_mg = b_rg_0:b_rg_1;
mag = 10*log10(Img(b_az,b_rg_0:b_rg_1));

%MagArr(3,9) = MagArr(3,9) - 3 ;
%MagArr(3,3) = MagArr(3,3) - 3 ;
FP07_LinePlot(2,idx_mg,mag,MagArr(3,8),MagArr(3,11),subtt,1,0, MagArr(3,7), MagArr(3,9), MagArr(3,15), MagArr(3,17), MagArr(3,18))

subtt = sprintf('$P_2$ PSLR = %.1f dB',MagArr(3,16));
title(subtt, 'Units', 'normalized', 'Position', [0.1, 0.1, 0],'horizontalAlignment', 'left','interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','k');

print(h_Fig, '-dpng','-r600','Figure11')
movefile('Figure11.png','Figures')


%Subplot
Scale = 1.2;
h_Fig=figure('PaperPositionMode', 'manual','PaperUnits','inches','PaperPosition',[0 0 3.5*2 3.5*2/1.618*Scale],'Position',[200 300 800 800/1.618*Scale]);
t = tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');

nexttile %subplot(3,2,1)
idx_mg =  1:etaTotal;
mag = 10*log10(Img(idx_mg,c_rg));

FP071_LinePlot(1,idx_mg,mag,MagArr(4,8),MagArr(4,11),subtt,0,0)
subtt = sprintf('$L_3$');
title(subtt, 'Units', 'normalized', 'Position', [0.5, 0.1, 0],'horizontalAlignment', 'left','interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','k');


nexttile %subplot(3,2,2)
idx_mg = c_az_0:c_az_1;
mag = 10*log10(Img(idx_mg,c_rg));

FP071_LinePlot(2,idx_mg,mag,MagArr(4,8),MagArr(4,11),subtt,0,0, MagArr(4,7), MagArr(4,9), MagArr(4,15), MagArr(4,17), MagArr(4,18))
%subtt = sprintf('$P_3$ PSLR = %.1f dB',MagArr(4,16));
subtt = sprintf('$P_3$');
title(subtt, 'Units', 'normalized', 'Position', [0.5, 0.1, 0],'horizontalAlignment', 'left','interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','k');


nexttile %subplot(3,2,3)
idx_mg = 1:etaTotal;
mag = 10*log10(Img(idx_mg,grp_rg));

FP07_LinePlot(1,idx_mg,mag,MagArr(5,8),MagArr(5,11),subtt,0,1)
subtt = sprintf('$LGRP_a$');
title(subtt, 'Units', 'normalized', 'Position', [0.45, 0.1, 0],'horizontalAlignment', 'left','interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','k');


nexttile %subplot(3,2,4)
idx_mg = grp_az_0:grp_az_1;
mag = 10*log10(Img(idx_mg,grp_rg));

FP071_LinePlot(2,idx_mg,mag,MagArr(5,8),MagArr(5,11),subtt,0,0, MagArr(5,7), MagArr(5,9), MagArr(5,15), MagArr(5,17), MagArr(5,18))
%subtt = sprintf('$GRP_{a}$ PSLR = %.1f dB',MagArr(5,16));
subtt = sprintf('$GRP_{a}$');
title(subtt, 'Units', 'normalized', 'Position', [0.45, 0.1, 0],'horizontalAlignment', 'left','interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','k');


nexttile %subplot(3,2,5)
idx_mg = 1:etaTotal;
mag = 10*log10(Img(idx_mg,d_rg));

FP071_LinePlot(1,idx_mg,mag,MagArr(6,8),MagArr(6,11),subtt,2,0)
subtt = sprintf('$L_4$');
title(subtt, 'Units', 'normalized', 'Position', [0.5, 0.1, 0],'horizontalAlignment', 'left','interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','k');

nexttile %subplot(3,2,6)
idx_mg = d_az_0:d_az_1;
mag = 10*log10(Img(idx_mg,d_rg));
%MagArr(6,9) = MagArr(6,9)-9;
%MagArr(6,3) = MagArr(6,3)-3;

FP071_LinePlot(2,idx_mg,mag,MagArr(6,8),MagArr(6,11),subtt,2,0, MagArr(6,7), MagArr(6,9), MagArr(6,15), MagArr(6,17), MagArr(6,18))
%subtt = sprintf('$P_4$ PSLR = = %.1f dB',MagArr(6,16));
subtt = sprintf('$P_4$');
title(subtt, 'Units', 'normalized', 'Position', [0.5, 0.1, 0],'horizontalAlignment', 'left','interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','k');

% %% left side
% idx_mg_lf = d_az_0:MagArr(6,2);
% mag_lf = 10*log10(Img(idx_mg_lf,d_az)); 
% [NuAmp, NuTime] = findpeaks(max(mag_lf)-mag_lf,idx_mg_lf);
% [~,idx] = sort(NuAmp,'ascend');
% MagArr(6,1) = NuTime(idx(1));%first null 
% MagArr(6,4) = max(mag)-NuAmp(idx(1));%left-side null 
% 
% %interpolation for plotting left-side null
% int_id = idx_mg_lf(1):0.001:idx_mg_lf(end);
% xSpline1 = interp1(idx_mg_lf,mag_lf,int_id,'spline');
% [NuAmp, NuTime] = findpeaks(max(xSpline1)-xSpline1,int_id);
% [~,idxp] = sort(NuAmp,'ascend');
% MagArr(6,7) = NuTime(idxp(1));%first null  
% MagArr(6,10) = max(xSpline1)-NuAmp(idxp(1));%left-side null 
% 
% %% right side
% idx_mg_rg = MagArr(6,2):a_az_1;
% mag_rg = 10*log10(Img(a_az,idx_mg_rg)); 
% [NuAmp, NuTime] = findpeaks(max(mag_rg)-mag_rg,idx_mg_rg);
% [~,idx] = sort(NuAmp,'ascend');
% MagArr(1,3) = NuTime(idx(1));%mainlobe 
% MagArr(1,6) = max(mag)-NuAmp(idx(1));%left-side null 
% 
% %interpolation for plotting right-side null
% int_id = idx_mg_rg(1):0.001:idx_mg_rg(end);
% xSpline1 = interp1(idx_mg_rg,mag_rg,int_id,'spline');
% [NuAmp, NuTime] = findpeaks(max(xSpline1)-xSpline1,int_id);
% [~,idxp] = sort(NuAmp,'ascend');
% MagArr(1,9) = NuTime(idxp(1));%mainlobe  
% MagArr(1,12) = max(xSpline1)-NuAmp(idxp(1));%left-side null 
% 
% hold on
% plot(MagArr(6,1:3),MagArr(6,4:6),'r*')

% %%
% %interpolation for plotting right-side null
% figure
% idx_mg = d_rg_0:d_rg_1;
% mag = 10*log10(Img(idx_mg,a_az));%%azimuth
% xSpline = interp1(idx_mg,mag,idx_mg(1):0.001:idx_mg(end),'spline');
% plot(idx_mg(1):0.001:idx_mg(end),xSpline,'-')
% hold on
% [NuAmp, NuTime] = findpeaks(max(xSpline1)-xSpline1,int_id);
% [~,idxp] = sort(NuAmp,'descend');
% hold on
% %plot(NuTime,NuAmp,'r*');
% %MagArr(1,9) = NuTime(idxp(1));%mainlobe  
% %MagArr(1,12) = max(xSpline1)-NuAmp(idxp(1));%left-side null 
% %%

%correction
%MagArr(4,1) = 

print(h_Fig, '-dpng','-r600','Figure10')
movefile('Figure10.png','Figures')

% figure
% s = surf(interp2(Img(hbin4y(6,3):hbin4y(6,4),hbin4x(6,3):hbin4x(6,4))));
% s.EdgeColor = 'none';
% xlabel('Range');
% ylabel('Azimuth');
% zlabel('Amplitude');
% %%%%%%%%%

