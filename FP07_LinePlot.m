function FP07_LinePlot(type, idx_mg,mag,o_ptsx,o_ptsy,subtt,xl_flag, yl_flag, shd0, shd1,p_sl,ymax, ymin)

if type == 1
xSpline = interp1(idx_mg,mag,idx_mg(1):0.001:idx_mg(end),'spline');
plot(idx_mg(1):0.001:idx_mg(end),xSpline,'-')
hold on
plot(o_ptsx,o_ptsy,'bO','MarkerSize',10)
xlim tight
ylim([-150 o_ptsy+2])
yticks(o_ptsy-30:10:o_ptsy)


if xl_flag==1
xlabel('Fast Time','interpreter','Tex')
elseif xl_flag==2
xlabel('Slow Time','interpreter','Tex')
else 
xticklabels('')
end


if yl_flag==1
ylabel('Magnitude[dB]','interpreter','Tex')
yticklabels(-30:10:0)
else 
yticklabels('')
end

%subtt = sprintf('(a)');
u = title(subtt, 'Units', 'normalized', 'Position', [0.5, 0.1, 0],'horizontalAlignment', 'left','interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','none');
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',10);
end 

if type == 2
    xSpline = interp1(idx_mg,mag,idx_mg(1):0.001:idx_mg(end),'spline');
    plot(idx_mg(1):0.001:idx_mg(end),xSpline,'-')
    hold on
    xlim tight

    ylim([-150 o_ptsy+2])
    yticks(o_ptsy-30:10:o_ptsy)
    yticklabels('')

    if xl_flag==1
        xlabel('Fast Time','interpreter','Tex')
    elseif xl_flag==2
        xlabel('Slow Time','interpreter','Tex')
    else 
        xticklabels('')
    end

    % shading
    x_points = [shd0, shd0, shd1, shd1];  
    y_points = [min(ylim), max(ylim), max(ylim), min(ylim)];
    color = [0, 0, 1];
    hold on;
    a = fill(x_points, y_points, color);
    a.FaceAlpha = 0.1;
    a.EdgeColor = "none";

    %average
    hold on;
    yline(p_sl,'-.r');
   % hold on
  %  yline([ymax ymin],'--');
   
title(subtt, 'Units', 'normalized', 'Position', [0.5, 0.1, 0],'horizontalAlignment', 'center','interpreter','latex','FontSize',10,'BackgroundColor','w','EdgeColor','k');
drawnow
set(gca,'LooseInset',get(gca,'TightInset'),'FontSize',10);

end
end