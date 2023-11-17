function[G,Tdp]=ConductiveHeatFluxFR_GroundImp(TemperatureC,TempDamp_ittm,TempVec_ittm,Owater_ittm,...
	ParCalculation,ParThermalGround,FractionsGround,ParSoilGround,ParVegTree,ParVegGround)

% TemperatureC(:,1)		=	Temperature ground impervious area
% TemperatureC(:,2)		=	Temperature ground bare area
% TemperatureC(:,3)		=	Temperature ground vegetated area
% TemperatureC(:,4)		=	Temperature sunlit area
% TemperatureC(:,5)		=	Temperature shaded area
% TemperatureC(:,6)		=	Temperature tree canopy
% TemperatureC(:,7)		=	Interior temperature sunlit wall
% TemperatureC(:,8)		=	Interior temperature shaded wall
% TemperatureC(:,9)		=	Temperature canyon
% TemperatureC(:,10)	=	specific humidity canyon

Ts				=	TemperatureC(1,1);	
Tdptm1			=	TempDamp_ittm.TDampGroundImp;
Tstm1			=	TempVec_ittm.TGroundImp;
dts				=	ParCalculation.dts;
lan_dry_imp		=	ParThermalGround.lan_dry_imp;
cv_s_imp		=	ParThermalGround.cv_s_imp;
Cimp			=	(FractionsGround.fimp>0);

% To calculate the influence of soil moisture on the soil heat capacity
Otm1		=	Owater_ittm.OwGroundSoilImp;
Pcla		=	ParSoilGround.Pcla;
Psan		=	ParSoilGround.Psan;
Porg		=	ParSoilGround.Porg;
Kfc			=	ParSoilGround.Kfc;
Phy			=	ParSoilGround.Phy;
SPAR		=	ParSoilGround.SPAR;
Kbot		=	ParSoilGround.Kbot;
CASE_ROOT_H	=	ParVegTree.CASE_ROOT;
CASE_ROOT_L	=	ParVegGround.CASE_ROOT;
ZR95_H		=	ParVegTree.ZR95;
ZR95_L		=	ParVegGround.ZR95;
ZR50_H		=	ParVegTree.ZR50;
ZR50_L		=	ParVegGround.ZR50;
ZRmax_H		=	ParVegTree.ZRmax;
ZRmax_L		=	ParVegGround.ZRmax;
Zs			=	ParSoilGround.Zs;

% Calculate soil properties underneath impervious surface
[~,dz,~,Osat,Ohy,~,~,~,~,~,~,~,~,~,...
	~,~,~,~,~,~,~,~,rsd,lan_dry,lan_s,cv_s]...
	=soil_functions.SoilParametersTotal(Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,Zs);

NotImp = ~isnan(Otm1) & Otm1~=0;

Otm1Ave = sum((dz(NotImp)./sum(dz(NotImp))).*Otm1(NotImp));

[lanS,cv_Soil,~]=soil_functions.Soil_Thermal_properties(Tdptm1-273.15,rsd,lan_dry,lan_s,cv_s,Osat,Ohy,Otm1Ave);

%--------------------------------------------------------------
lan_tot = (lan_dry_imp*sum(dz(~NotImp)) + lanS*sum(dz(NotImp)))/sum(dz);
cv_tot  = (cv_s_imp*sum(dz(~NotImp)) + cv_Soil*sum(dz(NotImp)))/sum(dz);   

tau		=	86400; %% [s] time constant
CTt		=	2*(sqrt(pi/(lan_tot*cv_tot*tau))); %%  [K m^2/J] Total Thermal Capacity Soil

[G,Tdp]	=	conductive_heat_functions.Soil_Heat(dts,Ts-273.15,Tstm1-273.15,Tdptm1-273.15,CTt);
Tdp		=	Tdp + 273.15;
G		=	G*Cimp;
