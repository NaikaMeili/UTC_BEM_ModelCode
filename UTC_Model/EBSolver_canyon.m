function[Ycanyon,G2WallSun,G2WallShade,SWRabs_t]=EBSolver_canyon(TemperatureC,TemperatureB,TempVec_ittm,Humidity_ittm,MeteoData,...
		Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,CiCO2Leaf_ittm,TempDamp_ittm,ViewFactor,...
		Gemeotry_m,ParTree,geometry,FractionsGround,...
		WallLayers,ParSoilGround,ParInterceptionTree,...
		PropOpticalGround,PropOpticalWall,PropOpticalTree,...
		ParThermalGround,ParThermalWall,ParVegGround,ParVegTree,...
		SunPosition,HumidityAtm,Anthropogenic,ParCalculation,...
        TempVecB_ittm,G2Roof,PropOpticalIndoors,ParHVAC,ParThermalBulidFloor,...
        ParWindows,BEM_on,RESPreCalc,fconvPreCalc,fconv,rsGroundPreCalc,rsTreePreCalc,HVACSchedule)


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
[SWRin_t,SWRout_t,SWRabs_t,SWRabsDir_t,SWRabsDiff_t,SWREB_t]...
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
[G1GroundImp,Tdp_ground_imp]=conductive_heat_functions.ConductiveHeatFluxFR_GroundImp(TemperatureC,TempDamp_ittm,TempVec_ittm,Owater_ittm,...
	ParCalculation,ParThermalGround,FractionsGround,ParSoilGround,ParVegTree,ParVegGround);

% Conductive heat flux bare ground
type	=	0; % bare
[G1GroundBare,Tdp_ground_bare]=conductive_heat_functions.ConductiveHeatFluxFR_GroundVegBare(TemperatureC,TempDamp_ittm,Owater_ittm,TempVec_ittm,...
	ParCalculation,ParSoilGround,ParVegGround,ParVegTree,FractionsGround,type);

% Conductive heat flux vegetated ground
type	=	1; % vegetated
[G1GroundVeg,Tdp_ground_veg]=conductive_heat_functions.ConductiveHeatFluxFR_GroundVegBare(TemperatureC,TempDamp_ittm,Owater_ittm,TempVec_ittm,...
	ParCalculation,ParSoilGround,ParVegGround,ParVegTree,FractionsGround,type);


%% Sensible and latent heat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[HfluxCanyon,LEfluxCanyon,ra_canyon,ra_orig,fconv,HumidityCan]=turbulent_heat_function.HeatFlux_canyon...
	(TemperatureC,Gemeotry_m,MeteoData,ParVegTree,ParTree,fconvPreCalc,fconv);

% Turbulent heat fluxes from ground and trees to canyon
[HfluxGroundImp,HfluxGroundBare,HfluxGroundVeg,HfluxTree,...
	Eground_imp_pond,Eground_bare_pond,Eground_bare_soil,Eground_veg_int,...
	Eground_veg_pond,Eground_veg_soil,TEground_veg,E_tree_int,TE_tree,...
	Ebare,Eveg,Etree,...
	LEfluxGroundImp,LEfluxGroundBarePond,LEfluxGroundBareSoil,LEfluxGroundVegInt,...
	LEfluxGroundVegPond,LEfluxGroundVegSoil,LTEfluxGroundVeg,LEfluxTreeInt,LTEfluxTree,...
	LEbare,LEveg,LEtree,...
	Ci_sun_tree,Ci_shd_tree,Ci_sun_ground,Ci_shd_ground,...
	rap_can,rap_Htree_In,rb_H,rb_L,...
	r_soil_bare,r_soil_veg,alp_soil_bare,alp_soil_veg,...
	rs_sun_L,rs_shd_L,rs_sun_H,rs_shd_H,...
	u_Hcan,u_Zref_und,Fsun_L,Fshd_L,dw_L]...
	=turbulent_heat_function.HeatFlux_ground(TemperatureC,TempVec_ittm,MeteoData,Gemeotry_m,geometry,FractionsGround,ParTree,...
	ParVegGround,ParVegTree,ParSoilGround,SoilPotW_ittm,Owater_ittm,Vwater_ittm,ExWater_ittm,Int_ittm,CiCO2Leaf_ittm,...
	ParInterceptionTree,ParCalculation,SWRabsDir_t.SWRabsTree,SWRabsDiff_t.SWRabsTree,SWRabsDir_t.SWRabsGroundVeg,SWRabsDiff_t.SWRabsGroundVeg,...
    RESPreCalc,rsGroundPreCalc,rsTreePreCalc);

% Turbulent heat fluxes from sunlit and shaded wall to canyon
[HfluxWallSun,HfluxWallShade,Ewsun,Ewshade,LEwsun,LEwshade,RES_w1,RES_w2,rap_Zp1_In,rap_Zp2_In,...
	Hwsun1,Hwshade1,Hwsun2,Hwshade2,cp_atm,rho_atm,L_heat,Zp1,Zp2,rap_Zp1,rap_Zp2]...
	=turbulent_heat_function.HeatFlux_wall(TemperatureC,Gemeotry_m,MeteoData,ParVegTree,ParTree,ParVegGround,FractionsGround);



% Anthropogenic heat input from building energy model
%--------------------------------------------------------------------------
if BEM_on==1
    SWRinWsun = SWRabs_t.SWRabsWallSunTransmitted;
    SWRinWshd = SWRabs_t.SWRabsWallShadeTransmitted;
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


%--------------------------------------------------------------------------
if ParTree.trees==0
SWRabs_t.SWRabsTree	=	0;
LWRabs_t.LWRabsTree	=	0;
end

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
                Ctree*4*geometry.radius_tree*HfluxTree - HfluxCanyon - dS_H_air ...
                + WasteHeat.SensibleFromVent_Can + WasteHeat.SensibleFromAC_Can + WasteHeat.SensibleFromHeat_Can;

Ycanyon(10)	=	Cimp*FractionsGround.fimp*LEfluxGroundImp + Cbare*FractionsGround.fbare*(LEfluxGroundBarePond+LEfluxGroundBareSoil) + ...
				Cveg*FractionsGround.fveg*(LEfluxGroundVegInt + LEfluxGroundVegPond + LEfluxGroundVegSoil + LTEfluxGroundVeg) +...
                Ctree*4*geometry.radius_tree*(LEfluxTreeInt + LTEfluxTree) - LEfluxCanyon - dS_LE_air...
                + WasteHeat.LatentFromVent_Can + WasteHeat.LatentFromAC_Can +  WasteHeat.LatentFromHeat_Can;

 



			
