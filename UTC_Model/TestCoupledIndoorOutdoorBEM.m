load('CalculationSGTest.mat')


%--------------------------------------------------------------------------
figure
tiledlayout(2,2)
nexttile
plot(MeteoDataRaw.Date(1:n),TempVecB.Tceiling-273.15,'b','Linewidth',1,'DisplayName','Ceiling')
hold on
plot(MeteoDataRaw.Date(1:n),TempVecB.Tinwallsun-273.15,'r','Linewidth',1,'DisplayName','Wall sun')
plot(MeteoDataRaw.Date(1:n),TempVecB.Tinwallshd-273.15,'m','Linewidth',1,'DisplayName','Wall shade')
plot(MeteoDataRaw.Date(1:n),TempVecB.Tinground-273.15,'g','Linewidth',1,'DisplayName','Ground')
plot(MeteoDataRaw.Date(1:n),TempVecB.Tbin-273.15,'k--','Linewidth',1,'DisplayName','Indoor air')
legend

nexttile
plot(MeteoDataRaw.Date(1:n),TempVecB.Tbin-273.15,'k','Linewidth',1,'DisplayName','Indoor air')
hold on
plot(MeteoDataRaw.Date(1:n),TempVec.TCanyon(1:n,1,1)-273.15,'b','Linewidth',1,'DisplayName','Outdoor air')
legend

nexttile
plot(MeteoDataRaw.Date(1:n),TempVec.TCanyon-273.15,'k','Linewidth',1,'DisplayName','Canyon')
hold on
plot(MeteoDataRaw.Date(1:n),TempVec.TWallSun-273.15,'r','Linewidth',1,'DisplayName','W_{sun,out}')
plot(MeteoDataRaw.Date(1:n),TempVec.TWallIntSun-273.15,'r--','Linewidth',1,'DisplayName','W_{sun,int}')
plot(MeteoDataRaw.Date(1:n),TempVecB.Tinwallsun-273.15,'m','Linewidth',1,'DisplayName','W_{sun,b}')
%
plot(MeteoDataRaw.Date(1:n),TempVec.TWallShade-273.15,'b','Linewidth',1,'DisplayName','W_{shd,out}')
plot(MeteoDataRaw.Date(1:n),TempVec.TWallIntShade-273.15,'b--','Linewidth',1,'DisplayName','W_{shd,int}')
plot(MeteoDataRaw.Date(1:n),TempVecB.Tinwallshd-273.15,'c','Linewidth',1,'DisplayName','W_{shd,b}')
legend

nexttile
plot(MeteoDataRaw.Date(1:n),TempVec.Tatm-273.15,'k','Linewidth',1,'DisplayName','T_{atm}')
hold on
plot(MeteoDataRaw.Date(1:n),TempVec.TRoofImp-273.15,'r','Linewidth',1,'DisplayName','R_{imp,out}')
plot(MeteoDataRaw.Date(1:n),TempVec.TRoofIntImp-273.15,'r--','Linewidth',1,'DisplayName','R_{imp,int}')
plot(MeteoDataRaw.Date(1:n),TempVecB.Tceiling-273.15,'m','Linewidth',1,'DisplayName','R_{b}')
legend

%--------------------------------------------------------------------------
figure
plot(MeteoDataRaw.Date,Gflux.G1WallSun(:,1,1),'r','Linewidth',1,'DisplayName','G1_{w,sun}')
hold on
plot(MeteoDataRaw.Date,Gflux.G2WallSun(:,1,1),'r--','Linewidth',1,'DisplayName','G2_{w,sun}')
plot(MeteoDataRaw.Date,Gflux.G1WallShade(:,1,1),'m','Linewidth',1,'DisplayName','G1_{w,shade}')
plot(MeteoDataRaw.Date,Gflux.G2WallShade(:,1,1),'m--','Linewidth',1,'DisplayName','G2_{w,shade}')
plot(MeteoDataRaw.Date,Gflux.G1Roof(:,1,1),'b','Linewidth',1,'DisplayName','G1_{roof}')
plot(MeteoDataRaw.Date,Gflux.G2Roof(:,1,1),'b--','Linewidth',1,'DisplayName','G2_{roof}')
legend




% figure
% tiledlayout(2,3);
% nexttile
% hold on
% % plot(MeteoDataRaw.Date(1:n),HbuildInt.HBinRoof,'DisplayName','Hb_{roof}')
% % plot(MeteoDataRaw.Date(1:n),HbuildInt.HbinWallSun,'DisplayName','Hb_{w,sun}')
% % plot(MeteoDataRaw.Date(1:n),HbuildInt.HbinWallshd,'DisplayName','Hb_{w,shade}')
% % plot(MeteoDataRaw.Date(1:n),HbuildInt.HBinGround,'DisplayName','Hb_{ground}')
% plot(MeteoDataRaw.Date(1:n),HbuildInt.HbuildInSurf,'DisplayName','Hb_{surftot}')
% plot(MeteoDataRaw.Date(1:n),HbuildInt.Hvent,'DisplayName','Hb_{ventilation}')
% plot(MeteoDataRaw.Date(1:n),HbuildInt.H_AC,'DisplayName','Hb_{AC}')
% plot(MeteoDataRaw.Date(1:n),HbuildInt.H_air,'DisplayName','Hb_{air,storage}')
% legend







% figure
% plot(MeteoDataRaw.Date(1:n),TempVecB.qbin,'DisplayName','Indoor air q')
% hold on
% plot(MeteoDataRaw.Date(1:n),Humidity.CanyonSpecific(1:n,1,1),'DisplayName','Outdoor air')
% legend
