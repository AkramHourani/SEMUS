function [pt_rg,pt_az,pt_rg_0,pt_rg_1,pt_az_0,pt_az_1] = FP09_PlotMag(rg_pt,az_pt,trg_val,max_rg,max_az)

iy = find(trg_val(:,1) == rg_pt & trg_val(:,2)== az_pt);
pt_rg = trg_val(iy,3);
pt_az = trg_val(iy,4);

iy = find(trg_val(:,1) == rg_pt-1 & trg_val(:,2)== az_pt);
pt_rg_0 = round(pt_rg-(pt_rg-trg_val(iy,3))/2);
if isempty(pt_rg_0) == 1 %%for edges, directly divide by 2
    pt_rg_0 = round(pt_rg/2);
end

iy = find(trg_val(:,1) == rg_pt+1 & trg_val(:,2)== az_pt);
pt_rg_1 = round(pt_rg+(trg_val(iy,3)-pt_rg)/2);
if isempty(pt_rg_1) == 1 %%for edges, directly divide by 2
    pt_rg_1 = round(max_rg-pt_rg/2);
end

iy = find(trg_val(:,1) == rg_pt & trg_val(:,2)== az_pt-1);
pt_az_0 = round(pt_az-(pt_az-trg_val(iy,4))/2);
if isempty(pt_az_0) == 1 %%for edges, directly divide by 2
    pt_az_0 = round(pt_az/2);
end

iy = find(trg_val(:,1)== rg_pt & trg_val(:,2)== az_pt+1);
pt_az_1 = round(pt_az+(trg_val(iy,4)-pt_az)/2);
if isempty(pt_az_1) == 1 %%for edges, directly divide by 2
    pt_az_1 = round(max_az-pt_az/2);
end

end