function[SWRin_t,SWRout_t,SWRabs_t,SWRabsDir_t,SWRabsDiff_t,SWREB_t,albedo_canyon,...
		LWRin_t,LWRout_t,LWRabs_t,LWREB_t,...
		HfluxGroundImp,HfluxGroundBare,HfluxGroundVeg,HfluxTree,HfluxGround,...
		EfluxGroundImp,EfluxGroundBarePond,EfluxGroundBareSoil,EfluxGroundVegInt,...
		EfluxGroundVegPond,EfluxGroundVegSoil,TEfluxGroundVeg,EfluxTreeInt,TEfluxTree,...
		EfluxGroundBare,EfluxGroundVeg,EfluxGround,EfluxTree,...
		LEfluxGroundImp,LEfluxGroundBarePond,LEfluxGroundBareSoil,LEfluxGroundVegInt,...
		LEfluxGroundVegPond,LEfluxGroundVegSoil,LTEfluxGroundVeg,LEfluxTreeInt,LTEfluxTree,...
		LEfluxGroundBare,LEfluxGroundVeg,LEfluxGround,LEfluxTree,...
		CiCO2LeafTreeSun,CiCO2LeafTreeShd,CiCO2LeafGroundVegSun,CiCO2LeafGroundVegShd,...
		raCanyontoAtm,raCanyontoAtmOrig,rap_can,rap_Htree_In,rb_HGround,rb_LGround,...
		r_soilGroundbare,r_soilGroundveg,alp_soil_bare,alp_soil_veg,...
		rs_sunGround,rs_shdGround,rs_sunTree,rs_shdTree,...
		Fsun_L,Fshd_L,dw_L,RES_w1,RES_w2,rap_W1_In,rap_W2_In,rap_Zp1,...
		HfluxWallSun,HfluxWallShade,EfluxWallSun,EfluxWallShade,LEfluxWallSun,LEfluxWallShade,HfluxCanyon,LEfluxCanyon,EfluxCanyon,...
		G1WallSun,G2WallSun,dsWallSun,G1WallShade,G2WallShade,dsWallShade,...
		G1GroundImp,TDampGroundImp,G1GroundBare,TDampGroundBare,G1GroundVeg,TDampGroundVeg,GTree,TDampTree,G1Ground,G1Canyon,G2Canyon,...
		dsGroundImp,dsGroundBare,dsGroundVeg,dsTree,dsCanyonAir,Ycanyon,...
		q_tree_dwn,In_tree,dIn_tree_dt,q_gveg_dwn,In_gveg,dIn_gveg_dt,...
		q_gimp_runoff,In_gimp,dIn_gimp_dt,f_inf_gimp,q_gbare_runoff,In_gbare,dIn_gbare_dt,f_inf_gbare,...
		q_gveg_runoff,In_gveg_pond,dIn_gveg_pond_dt,f_inf_gveg,...
		V_gimp,O_gimp,OS_gimp,Lk_gimp,Psi_s_H_gimp,Psi_s_L_gimp,...
		Exwat_H_gimp,Exwat_L_gimp,Rd_gimp,TEgveg_imp,TEtree_imp,...
		Egimp_soil,dV_dt_gimp,Psi_Soil_gimp,Kf_gimp,...
		V_gbare,O_gbare,OS_gbare,Lk_gbare,Psi_s_H_gbare,Psi_s_L_gbare,...
		Exwat_H_gbare,Exwat_L_gbare,Rd_gbare,TEgveg_bare,TEtree_bare,...
		Egbare_Soil,dV_dt_gbare,Psi_soil_gbare,Kf_gbare,...
		V_gveg,O_gveg,OS_gveg,Lk_gveg,Psi_s_H_gveg,Psi_s_L_gveg,...
		Exwat_H_gveg,Exwat_L_gveg,Rd_gveg,TEgveg_veg,TEtree_veg,...
		Egveg_Soil,dV_dt_gveg,Psi_soil_gveg,Kf_gveg,...
		Qin_imp,Qin_bare,Qin_veg,Qin_bare2imp,Qin_bare2veg,Qin_imp2bare,Qin_imp2veg,Qin_veg2imp,Qin_veg2bare,...
		V,O,OS,Lk,Rd,dV_dt,Psi_s_L,Exwat_L,TEgveg_tot,Psi_s_H_tot,Exwat_H,...
		TEtree_tot,EB_TEtree,EB_TEgveg,WBIndv,WBTot,...
		Runoff,Runon_ittm,Etot,DeepGLk,StorageTot,...
		EBGroundImp,EBGroundBare,EBGroundVeg,EBTree,EBWallSun,EBWallShade,EBWallSunInt,EBWallShadeInt,EBCanyonT,EBCanyonQ,...
		HumidityCan,HumidityAtm,u_Hcan,u_Zref_und,T2m,q2m,e_T2m,RH_T2m,qcan,e_Tcan,RH_Tcan,...
		DHi,Himp_2m,Hbare_2m,Hveg_2m,Hwsun_2m,Hwshade_2m,Hcan_2m,...
		DEi,Eimp_2m,Ebare_soil_2m,Eveg_int_2m,Eveg_soil_2m,TEveg_2m,Ecan_2m,dS_H_air,dS_LE_air]...
		=EB_WB_canyon(TemperatureC,TemperatureB,TempVec_ittm,Humidity_ittm,MeteoData,...
		Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,CiCO2Leaf_ittm,...
        TempDamp_ittm,Runon_ittm,Qinlat_ittm,ViewFactor,...
		Gemeotry_m,ParTree,geometry,FractionsGround,...
		WallLayers,ParSoilGround,ParInterceptionTree,...
		PropOpticalGround,PropOpticalWall,PropOpticalTree,...
		ParThermalGround,ParThermalWall,ParVegGround,ParVegTree,...
		SunPosition,HumidityAtm,Anthropogenic,ParCalculation,...
        TempVecB_ittm,G2Roof,PropOpticalIndoors,ParHVAC,ParThermalBulidFloor,ParWindows,BEM_on,...
        RESPreCalc,fconvPreCalc,fconv,rsGroundPreCalc,rsTreePreCalc,HVACSchedule)

	
% Temperature vector:
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


%% Shortwave radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SWRin_t,SWRout_t,SWRabs_t,SWRabsDir_t,SWRabsDiff_t,SWREB_t,albedo_canyon]...
         =radiation_functions.TotalSWRabsorbed(geometry,FractionsGround,ParTree,...
		 PropOpticalGround,PropOpticalWall,PropOpticalTree,ParVegTree,MeteoData,...
         SunPosition,ViewFactor,ParWindows,BEM_on);
	 
% Tree absorbed: conversion from sphere to horizontal projected area
SWRabs_t.SWRabsTree		=	SWRabs_t.SWRabsTree*4*geometry.radius_tree*pi/(4*geometry.radius_tree);
SWRabsDir_t.SWRabsTree	=	SWRabsDir_t.SWRabsTree*4*geometry.radius_tree*pi/(4*geometry.radius_tree);
SWRabsDiff_t.SWRabsTree	=	SWRabsDiff_t.SWRabsTree*4*geometry.radius_tree*pi/(4*geometry.radius_tree);

if BEM_on==1
    if ParWindows.WindowsOn==0
        ParWindows.GlazingRatio = 0;
    end

    SWRabs_t.SWRabsWallSun      =   SWRabs_t.SWRabsWallSun; % W/m^2 wall
    SWRabs_t.SWRabsWindowSun    =   0; % W/m^2 window
    SWRabs_t.SWRtransWindowSun  =   SWRabs_t.SWRabsWallSun; % W/m^2 window
    SWRabs_t.SWRabsWallSunExt           =   (1-ParWindows.GlazingRatio).*SWRabs_t.SWRabsWallSun + ParWindows.GlazingRatio.*SWRabs_t.SWRabsWindowSun;
    SWRabs_t.SWRabsWallSunTransmitted   =   ParWindows.GlazingRatio.*SWRabs_t.SWRtransWindowSun;
    
    SWRabs_t.SWRabsWallShade      =   SWRabs_t.SWRabsWallShade; % W/m^2 wall
    SWRabs_t.SWRabsWindowShade    =   0; % W/m^2 window
    SWRabs_t.SWRtransWindowShade  =   SWRabs_t.SWRabsWallShade; % W/m^2 window
    SWRabs_t.SWRabsWallShadeExt           =   (1-ParWindows.GlazingRatio).*SWRabs_t.SWRabsWallShade + ParWindows.GlazingRatio.*SWRabs_t.SWRabsWindowShade;
    SWRabs_t.SWRabsWallShadeTransmitted   =   ParWindows.GlazingRatio.*SWRabs_t.SWRtransWindowShade;
else
    SWRabs_t.SWRabsWallSun      =   SWRabs_t.SWRabsWallSun; % W/m^2 wall
    SWRabs_t.SWRabsWindowSun    =   0; % W/m^2 window
    SWRabs_t.SWRtransWindowSun  =   0; % W/m^2 window
    SWRabs_t.SWRabsWallSunExt           =   SWRabs_t.SWRabsWallSun;
    SWRabs_t.SWRabsWallSunTransmitted   =   0;
    
    SWRabs_t.SWRabsWallShade      =   SWRabs_t.SWRabsWallShade; % W/m^2 wall
    SWRabs_t.SWRabsWindowShade    =   0; % W/m^2 window
    SWRabs_t.SWRtransWindowShade  =   0; % W/m^2 window
    SWRabs_t.SWRabsWallShadeExt   =   SWRabs_t.SWRabsWallShade;
    SWRabs_t.SWRabsWallShadeTransmitted   =   0;
end


if abs(SWRabs_t.SWRabsWallSunExt+SWRabs_t.SWRabsWallSunTransmitted - SWRabs_t.SWRabsWallSun) > 10^-8
    disp('Warning')
elseif abs(SWRabs_t.SWRabsWallShadeExt + SWRabs_t.SWRabsWallShadeTransmitted - SWRabs_t.SWRabsWallShade) > 10^-8 
    disp('Warning')
end


%% Longwave radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[LWRin_t,LWRout_t,LWRabs_t,LWREB_t]...
		 =radiation_functions.TotalLWRabsorbed(TemperatureC,geometry,MeteoData,...
		 FractionsGround,PropOpticalGround,PropOpticalWall,PropOpticalTree,ParTree,ViewFactor);

% Tree absorbed: conversion from sphere to horizontal projected area
LWRabs_t.LWRabsTree		=	LWRabs_t.LWRabsTree*4*geometry.radius_tree*pi/(4*geometry.radius_tree);

%% Conductive heat fluxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conductive heat flux sunlit wall
type	=	1;
[G1WallSun,G2WallSun,dsWallSun]=conductive_heat_functions.ConductiveHeatFlux_Walls...
	(TemperatureC,TemperatureB,TempVec_ittm,TempVecB_ittm,Anthropogenic,ParThermalWall,WallLayers,ParCalculation,type,ParWindows,BEM_on);

% Conductive heat flux shaded wall
type	=	0;
[G1WallShade,G2WallShade,dsWallShade]=conductive_heat_functions.ConductiveHeatFlux_Walls...
	(TemperatureC,TemperatureB,TempVec_ittm,TempVecB_ittm,Anthropogenic,ParThermalWall,WallLayers,ParCalculation,type,ParWindows,BEM_on);

% Conductive heat flux impervious ground
[G1GroundImp,TDampGroundImp]=conductive_heat_functions.ConductiveHeatFluxFR_GroundImp(TemperatureC,TempDamp_ittm,TempVec_ittm,Owater_ittm,...
	ParCalculation,ParThermalGround,FractionsGround,ParSoilGround,ParVegTree,ParVegGround);

% Conductive heat flux bare ground
type	=	0; % bare
[G1GroundBare,TDampGroundBare]=conductive_heat_functions.ConductiveHeatFluxFR_GroundVegBare(TemperatureC,TempDamp_ittm,Owater_ittm,TempVec_ittm,...
	ParCalculation,ParSoilGround,ParVegGround,ParVegTree,FractionsGround,type);

% Conductive heat flux vegetated ground
type	=	1; % vegetated
[G1GroundVeg,TDampGroundVeg]=conductive_heat_functions.ConductiveHeatFluxFR_GroundVegBare(TemperatureC,TempDamp_ittm,Owater_ittm,TempVec_ittm,...
	ParCalculation,ParSoilGround,ParVegGround,ParVegTree,FractionsGround,type);


G1Ground	=	G1GroundImp*FractionsGround.fimp + G1GroundBare*FractionsGround.fbare + G1GroundVeg*FractionsGround.fveg;
GTree		=	NaN;
dsGroundImp	=	NaN;
dsGroundBare=	NaN;
dsGroundVeg	=	NaN;
dsTree		=	NaN;
dsCanyonAir	=	NaN;
TDampTree	=	NaN;
G1Canyon	=	geometry.wcanyon*G1Ground + geometry.hcanyon*G1WallSun + geometry.hcanyon*G1WallShade;
G2Canyon	=	geometry.wcanyon*0 + geometry.hcanyon*G2WallSun + geometry.hcanyon*G2WallShade;

%% Sensible and latent heat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[HfluxCanyon,LEfluxCanyon,raCanyontoAtm,raCanyontoAtmOrig,fconv,HumidityCan]=turbulent_heat_function.HeatFlux_canyon...
	(TemperatureC,Gemeotry_m,MeteoData,ParVegTree,ParTree,fconvPreCalc,fconv);

% Turbulent heat fluxes from ground and trees to canyon
[HfluxGroundImp,HfluxGroundBare,HfluxGroundVeg,HfluxTree,...
	EfluxGroundImp,EfluxGroundBarePond,EfluxGroundBareSoil,EfluxGroundVegInt,...
	EfluxGroundVegPond,EfluxGroundVegSoil,TEfluxGroundVeg,EfluxTreeInt,TEfluxTree,...
	Ebare,Eveg,Etree,...
	LEfluxGroundImp,LEfluxGroundBarePond,LEfluxGroundBareSoil,LEfluxGroundVegInt,...
	LEfluxGroundVegPond,LEfluxGroundVegSoil,LTEfluxGroundVeg,LEfluxTreeInt,LTEfluxTree,...
	LEbare,LEveg,LEtree,...
	CiCO2LeafTreeSun,CiCO2LeafTreeShd,CiCO2LeafGroundVegSun,CiCO2LeafGroundVegShd,...
	rap_can,rap_Htree_In,rb_HGround,rb_LGround,...
	r_soilGroundbare,r_soilGroundveg,alp_soil_bare,alp_soil_veg,...
	rs_sunGround,rs_shdGround,rs_sunTree,rs_shdTree,...
	u_Hcan,u_Zref_und,Fsun_L,Fshd_L,dw_L]...
	=turbulent_heat_function.HeatFlux_ground(TemperatureC,TempVec_ittm,MeteoData,Gemeotry_m,geometry,FractionsGround,ParTree,...
	ParVegGround,ParVegTree,ParSoilGround,SoilPotW_ittm,Owater_ittm,Vwater_ittm,ExWater_ittm,Int_ittm,CiCO2Leaf_ittm,...
	ParInterceptionTree,ParCalculation,SWRabsDir_t.SWRabsTree,SWRabsDiff_t.SWRabsTree,SWRabsDir_t.SWRabsGroundVeg,SWRabsDiff_t.SWRabsGroundVeg,...
    RESPreCalc,rsGroundPreCalc,rsTreePreCalc);

% Turbulent heat fluxes from sunlit and shaded wall to canyon
[HfluxWallSun,HfluxWallShade,EfluxWallSun,EfluxWallShade,LEfluxWallSun,LEfluxWallShade,RES_w1,RES_w2,rap_W1_In,rap_W2_In,...
	Hwsun1,Hwshade1,Hwsun2,Hwshade2,cp_atm,rho_atm,L_heat,Zp1,Zp2,rap_Zp1,rap_Zp2]...
	=turbulent_heat_function.HeatFlux_wall(TemperatureC,Gemeotry_m,MeteoData,ParVegTree,ParTree,ParVegGround,FractionsGround);



HfluxGround		=	HfluxGroundImp*FractionsGround.fimp + HfluxGroundBare*FractionsGround.fbare + HfluxGroundVeg*FractionsGround.fveg;

LEfluxGroundBare=	LEfluxGroundBarePond + LEfluxGroundBareSoil;
LEfluxGroundVeg	=	LEfluxGroundVegInt + LEfluxGroundVegPond + LEfluxGroundVegSoil + LTEfluxGroundVeg;
LEfluxGround	=	LEfluxGroundImp*FractionsGround.fimp + LEfluxGroundBare*FractionsGround.fbare + LEfluxGroundVeg*FractionsGround.fveg;
LEfluxTree		=	LEfluxTreeInt + LTEfluxTree;

EfluxGroundBare	=	EfluxGroundBarePond + EfluxGroundBareSoil;
EfluxGroundVeg	=	EfluxGroundVegInt + EfluxGroundVegPond + EfluxGroundVegSoil + TEfluxGroundVeg;
EfluxGround		=	EfluxGroundImp*FractionsGround.fimp + EfluxGroundBare*FractionsGround.fbare + EfluxGroundVeg*FractionsGround.fveg;
EfluxTree		=	EfluxTreeInt + TEfluxTree;
EfluxCanyon		=	NaN;

if BEM_on==1
    SWRinWsun = SWRabs_t.SWRabsWallSunTransmitted;
    SWRinWshd = SWRabs_t.SWRabsWallShadeTransmitted;
    
    % Building energy model to calculate waste heat into the canyon
    [~,WasteHeat]=BuildingEnergyModel.EBSolver_Building(TemperatureC,TemperatureB,TempVecB_ittm,TempVec_ittm,Humidity_ittm,MeteoData,...
        SWRinWsun,SWRinWshd,G2Roof,G2WallSun,G2WallShade,TempDamp_ittm,SWRabs_t,...
        Gemeotry_m,PropOpticalIndoors,ParHVAC,ParCalculation,ParThermalBulidFloor,ParWindows,BEM_on,HVACSchedule);
else
    WasteHeat.SensibleFromAC_Can = 0; WasteHeat.LatentFromAC_Can = 0;
    WasteHeat.WaterFromAC_Can = 0; WasteHeat.SensibleFromHeat_Can = 0; 
    WasteHeat.LatentFromHeat_Can = 0; WasteHeat.SensibleFromVent_Can = 0;
    WasteHeat.LatentFromVent_Can = 0;  WasteHeat.TotAnthInput_URB = 0;
end

% Change in heat storage in canyon air
%--------------------------------------------------------------------------
% TemperatureC(:,9)		=	Temperature canyon
% TemperatureC(:,10)	=	specific humidity canyon
Vcanyon = (Gemeotry_m.Width_canyon*Gemeotry_m.Height_canyon)/Gemeotry_m.Width_canyon;
dS_H_air = Vcanyon*cp_atm*rho_atm*(TemperatureC(:,9)-TempVec_ittm.TCanyon)/ParCalculation.dts; % (W/m), m^2*J/(kg K)*(kg/m^3)*K/s
dS_LE_air = Vcanyon*rho_atm*L_heat*(TemperatureC(:,10)-Humidity_ittm.CanyonSpecific)/ParCalculation.dts; % (W/m), m^2*(kg/m^3)*(J/kg)*(kg/kg)/s


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cimp			=	(FractionsGround.fimp>0);
Cbare			=	(FractionsGround.fbare>0);
Cveg			=	(FractionsGround.fveg>0);
Ctree			=	(ParTree.trees==1);

if FractionsGround.fimp>0
	Ycanyon(1)	=	SWRabs_t.SWRabsGroundImp + LWRabs_t.LWRabsGroundImp - G1GroundImp - HfluxGroundImp - LEfluxGroundImp;
else
	Ycanyon(1)	=	TemperatureC(1,1)-273.15;
end
if FractionsGround.fbare>0
	Ycanyon(2)	=	SWRabs_t.SWRabsGroundBare + LWRabs_t.LWRabsGroundBare - G1GroundBare - HfluxGroundBare - LEfluxGroundBarePond - LEfluxGroundBareSoil;
else
	Ycanyon(2)	=	TemperatureC(1,2)-273.15;
end
if FractionsGround.fveg>0
	Ycanyon(3)	=	SWRabs_t.SWRabsGroundVeg + LWRabs_t.LWRabsGroundVeg - G1GroundVeg - HfluxGroundVeg - LEfluxGroundVegInt - LEfluxGroundVegPond - LEfluxGroundVegSoil - LTEfluxGroundVeg;
else
	Ycanyon(3)	=	TemperatureC(1,3)-273.15;
end
Ycanyon(4)	=	SWRabs_t.SWRabsWallSunExt + LWRabs_t.LWRabsWallSun - G1WallSun - HfluxWallSun;
Ycanyon(5)	=	SWRabs_t.SWRabsWallShadeExt + LWRabs_t.LWRabsWallShade - G1WallShade - HfluxWallShade;
if ParTree.trees>0
	Ycanyon(6)	=	SWRabs_t.SWRabsTree + LWRabs_t.LWRabsTree - HfluxTree- LEfluxTreeInt - LTEfluxTree;		
else
	Ycanyon(6)	=	TemperatureC(1,6)-273.15;
end
Ycanyon(7)	=	G1WallSun - G2WallSun - dsWallSun; % Energy budget interior of sunlit wall
Ycanyon(8)	=	G1WallShade - G2WallShade - dsWallShade; % Energy budget interior of shaded wall
Ycanyon(9)	=	Anthropogenic.Qf_canyon + Cimp*FractionsGround.fimp*HfluxGroundImp + Cbare*FractionsGround.fbare*HfluxGroundBare + ...
                Cveg*FractionsGround.fveg*HfluxGroundVeg + geometry.hcanyon*HfluxWallSun + geometry.hcanyon*HfluxWallShade + ...
                Ctree*4*geometry.radius_tree*HfluxTree - HfluxCanyon - dS_H_air...
                + WasteHeat.SensibleFromVent_Can + WasteHeat.SensibleFromAC_Can + WasteHeat.SensibleFromHeat_Can;

Ycanyon(10)	=	Cimp*FractionsGround.fimp*LEfluxGroundImp + Cbare*FractionsGround.fbare*(LEfluxGroundBarePond+LEfluxGroundBareSoil) + ...
				Cveg*FractionsGround.fveg*(LEfluxGroundVegInt + LEfluxGroundVegPond + LEfluxGroundVegSoil + LTEfluxGroundVeg) +...
                Ctree*4*geometry.radius_tree*(LEfluxTreeInt + LTEfluxTree) - LEfluxCanyon - dS_LE_air...
                + WasteHeat.LatentFromVent_Can + WasteHeat.LatentFromAC_Can + WasteHeat.LatentFromHeat_Can;


if 	max(abs(Ycanyon(:)))>1
end
			
EBGroundImp		=	Ycanyon(1);
EBGroundBare	=	Ycanyon(2);
EBGroundVeg		=	Ycanyon(3);
EBTree			=	Ycanyon(6);
EBWallSun		=	Ycanyon(4);
EBWallShade		=	Ycanyon(5);
EBWallSunInt	=	Ycanyon(7);
EBWallShadeInt	=	Ycanyon(8);
EBCanyonT		=	Ycanyon(9);
EBCanyonQ		=	Ycanyon(10);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Air temperature at 2 m
[T2m,DHi,Himp_2m,Hbare_2m,Hveg_2m,Hwsun_2m,Hwshade_2m,Hcan_2m]=...
    turbulent_heat_function.CalculateT2m(TemperatureC(1,1),TemperatureC(1,2),TemperatureC(1,3),TemperatureC(1,4),TemperatureC(1,5),TemperatureC(1,9),...
	Zp1,rap_Zp1,rap_W1_In,rb_LGround,RES_w1,FractionsGround,Gemeotry_m,geometry,ParVegGround,...
    TempVec_ittm,cp_atm,rho_atm,ParCalculation,fconv,MeteoData);

% Air humidity at 2m
q2m	=fzero(@turbulent_heat_function.AirHumidity2m,TemperatureC(1,10),optimset(@fzero),T2m,...
	TemperatureC(1,1),TemperatureC(1,2),TemperatureC(1,3),TemperatureC(1,9),TemperatureC(1,10),...
	rap_Zp1,rap_W1_In,rb_LGround,alp_soil_bare,r_soilGroundbare,alp_soil_veg,r_soilGroundveg,rs_sunGround,rs_shdGround,...
	dw_L,Fsun_L,Fshd_L,FractionsGround,ParVegGround,...
	EfluxGroundImp,Ebare,EfluxGroundVegInt,EfluxGroundVegPond,EfluxGroundVegSoil,TEfluxGroundVeg,MeteoData.Pre,...
    Humidity_ittm,fconv,MeteoData,Gemeotry_m,rho_atm,Zp1,ParCalculation);

% Latent flux outputs
[DEi,Eimp_2m,Ebare_soil_2m,Eveg_int_2m,Eveg_soil_2m,TEveg_2m,Ecan_2m,...
	q2m,e_T2m,RH_T2m,qcan,e_Tcan,RH_Tcan...
	]=turbulent_heat_function.AirHumidity2mOutput(q2m,T2m,...
	TemperatureC(1,1),TemperatureC(1,2),TemperatureC(1,3),TemperatureC(1,9),TemperatureC(1,10),...
	rap_Zp1,rap_W1_In,rb_LGround,alp_soil_bare,r_soilGroundbare,alp_soil_veg,r_soilGroundveg,rs_sunGround,rs_shdGround,...
	dw_L,Fsun_L,Fshd_L,FractionsGround,ParVegGround,...
	EfluxGroundImp,Ebare,EfluxGroundVegInt,EfluxGroundVegPond,EfluxGroundVegSoil,TEfluxGroundVeg,MeteoData.Pre,...
    Humidity_ittm,fconv,MeteoData,Gemeotry_m,rho_atm,Zp1,ParCalculation);


%% Water fluxes and water balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[q_tree_dwn,In_tree,dIn_tree_dt,q_gveg_dwn,In_gveg,dIn_gveg_dt,...
	q_gimp_runoff,In_gimp,dIn_gimp_dt,f_inf_gimp,q_gbare_runoff,In_gbare,dIn_gbare_dt,f_inf_gbare,...
	q_gveg_runoff,In_gveg_pond,dIn_gveg_pond_dt,f_inf_gveg,...
	V_gimp,O_gimp,OS_gimp,Lk_gimp,Psi_s_H_gimp,Psi_s_L_gimp,...
	Exwat_H_gimp,Exwat_L_gimp,Rd_gimp,TEgveg_imp,TEtree_imp,...
	Egimp_soil,dV_dt_gimp,Psi_Soil_gimp,Kf_gimp,...
	V_gbare,O_gbare,OS_gbare,Lk_gbare,Psi_s_H_gbare,Psi_s_L_gbare,...
	Exwat_H_gbare,Exwat_L_gbare,Rd_gbare,TEgveg_bare,TEtree_bare,...
	Egbare_Soil,dV_dt_gbare,Psi_soil_gbare,Kf_gbare,...
	V_gveg,O_gveg,OS_gveg,Lk_gveg,Psi_s_H_gveg,Psi_s_L_gveg,...
	Exwat_H_gveg,Exwat_L_gveg,Rd_gveg,TEgveg_veg,TEtree_veg,...
	Egveg_Soil,dV_dt_gveg,Psi_soil_gveg,Kf_gveg,...
	Qin_imp,Qin_bare,Qin_veg,Qin_bare2imp,Qin_bare2veg,Qin_imp2bare,Qin_imp2veg,Qin_veg2imp,Qin_veg2bare,...
	V,O,OS,Lk,Rd,dV_dt,Psi_s_L,Exwat_L,TEgveg_tot,Psi_s_H_tot,Exwat_H,...
	TEtree_tot,EB_TEtree,EB_TEgveg,WBIndv,WBTot,...
	Runoff,Runon_ittm,Etot,DeepGLk,StorageTot]=...
	water_functions.WaterCanyon(MeteoData,Int_ittm,Owater_ittm,Runon_ittm,Qinlat_ittm,...
	EfluxTreeInt,EfluxGroundVegInt,EfluxGroundImp,EfluxGroundBarePond,EfluxGroundVegPond,EfluxGroundBareSoil,EfluxGroundVegSoil,TEfluxGroundVeg,TEfluxTree,...
	ParSoilGround,ParInterceptionTree,ParCalculation,ParVegGround,ParVegTree,...
	FractionsGround,geometry,ParTree,Gemeotry_m,Anthropogenic);

