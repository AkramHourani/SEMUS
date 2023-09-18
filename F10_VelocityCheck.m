%% This script to find the maximum Doppler frequency using the four corners of the swath 
function [V_max] = F10_VelocityCheck(latSwathL1SoI,lonSwathL1SoI,latSwathL2SoI,lonSwathL2SoI,SatllaSoI,E,Ro,Param)
[~,~,Rcr1] = geodetic2aer(latSwathL1SoI(1),lonSwathL1SoI(1),0,SatllaSoI(:,1),SatllaSoI(:,2),SatllaSoI(:,3),E);              % Reference range at the first corner of the near edge
[~,~,Rcr2] = geodetic2aer(latSwathL1SoI(end),lonSwathL1SoI(end),0,SatllaSoI(:,1),SatllaSoI(:,2),SatllaSoI(:,3),E);          % Reference range at the second corner of the near edge
[~,~,Rcr3] = geodetic2aer(latSwathL2SoI(1),lonSwathL2SoI(1),0,SatllaSoI(:,1),SatllaSoI(:,2),SatllaSoI(:,3),E);              % Reference range at the first corner of the far edge
[~,~,Rcr4] = geodetic2aer(latSwathL2SoI(end),lonSwathL2SoI(end),0,SatllaSoI(:,1),SatllaSoI(:,2),SatllaSoI(:,3),E);          % Reference range at the second corner of the far edge
Vo = diff(Ro)/Param.ts;                  % Velocity of the GRP
V1 = diff(Rcr1)/Param.ts;               % Velocity of the first corner of the near edge
V2 = diff(Rcr2)/Param.ts;               % Velocity of the second corner of the near edge
V3 = diff(Rcr3)/Param.ts;               % Velocity of the first corner of the far edge
V4 = diff(Rcr4)/Param.ts;               % Velocity of the second corner of the far edge
V = cat(1,Vo(:),V1(:),V2(:),V3(:),V4(:));
V_max = max(V,[],1);
% plot(abs([Vo,V1,V2,V3,V4]))
% PRFmin = 2 * RadPar.fo * V_max/ c
end