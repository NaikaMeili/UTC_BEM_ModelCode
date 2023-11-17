function[G,Tdp]=ConductiveHeatFluxFR_GroundVegBare(TemperatureC,TempDamp_ittm,Owater_ittm,TempVec_ittm,...
	ParCalculation,ParSoilGround,ParVegGround,ParVegTree,FractionsGround,type)

if type == 0 % bare soil
	Ts			=	TemperatureC(1,2);
	Tdptm1		=	TempDamp_ittm.TDampGroundBare;
	Otm1		=	Owater_ittm.OwGroundSoilBare;
	Tstm1		=	TempVec_ittm.TGroundBare;
	dts			=	ParCalculation.dts;
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
	Csoil		=	(FractionsGround.fbare>0);
elseif type ==1 % vegetated ground
	Ts			=	TemperatureC(1,3);
	Tdptm1		=	TempDamp_ittm.TDampGroundVeg;
	Otm1		=	Owater_ittm.OwGroundSoilVeg;
	Tstm1		=	TempVec_ittm.TGroundVeg;
	dts			=	ParCalculation.dts;
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
	Csoil		=	(FractionsGround.fveg>0);
else
	disp('please enter a valid specification of ground type. 0 = bare, 1 = vegetated')
end


[~,dz,~,Osat,Ohy,~,~,~,~,~,~,~,~,~,...
	~,~,~,~,~,~,~,~,rsd,lan_dry,lan_s,cv_s]...
	=soil_functions.SoilParametersTotal(Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,Zs);

Otm1Ave = sum((dz./sum(dz)).*Otm1);

[~,~,CTt]=soil_functions.Soil_Thermal_properties(Tdptm1-273.15,rsd,lan_dry,lan_s,cv_s,Osat,Ohy,Otm1Ave);

[G,Tdp]=conductive_heat_functions.Soil_Heat(dts,Ts-273.15,Tstm1-273.15,Tdptm1-273.15,CTt);
Tdp		=	Tdp + 273.15;
G		=	G*Csoil;
