%% This script to find the maximum Doppler frequency using the four corners of the swath 
function [V_max] = F03_VelocityCheck(latSwathL1,lonSwathL1,latSwathL2,lonSwathL2,Satlla,E,R,Param)
[~,~,Rcr1] = geodetic2aer(latSwathL1(1),lonSwathL1(1),0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);              % Reference range at the first corner of the near edge
[~,~,Rcr2] = geodetic2aer(latSwathL1(end),lonSwathL1(end),0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);          % Reference range at the second corner of the near edge
[~,~,Rcr3] = geodetic2aer(latSwathL2(1),lonSwathL2(1),0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);              % Reference range at the first corner of the far edge
[~,~,Rcr4] = geodetic2aer(latSwathL2(end),lonSwathL2(end),0,Satlla(:,1),Satlla(:,2),Satlla(:,3),E);          % Reference range at the second corner of the far edge
Vo = diff(R)/Param.tg;                  % Velocity of the GRP
V1 = diff(Rcr1)/Param.tg;               % Velocity of the first corner of the near edge
V2 = diff(Rcr2)/Param.tg;               % Velocity of the second corner of the near edge
V3 = diff(Rcr3)/Param.tg;               % Velocity of the first corner of the far edge
V4 = diff(Rcr4)/Param.tg;               % Velocity of the second corner of the far edge
V = cat(1,Vo(:),V1(:),V2(:),V3(:),V4(:));
V_max = max(V,[],1);
% plot(abs([Vo,V1,V2,V3,V4]))
% f_D_max = 2 * 2 * pi * V_max * RadPar.fo /c;
% f_D_max = -2 * 2 * pi * max(Vo) * RadPar.fo /c;
end
