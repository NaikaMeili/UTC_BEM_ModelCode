function[Ycanyon]=EBSolver_canyon(TemperatureC,TempVec,MeteoData,itt,...
		Int,ExWater,Vwater,Owater,SoilPotW,CiCO2Leaf,TempDamp,ViewFactor,...
		Gemeotry_m,ParTree,geometry,FractionsGround,...
		WallLayers,ParSoilGround,ParInterceptionTree,...
		PropOpticalGround,PropOpticalWall,PropOpticalTree,...
		ParThermalGround,ParThermalWall,ParVegGround,ParVegTree,...
		SunPosition,HumidityAtm,Anthropogenic,ParCalculation)

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
         SunPosition,ViewFactor);
	 
% Tree absorbed: conversion from sphere to horizontal projected area
SWRabs_t.SWRabsTree		=	SWRabs_t.SWRabsTree*4*geometry.radius_tree*pi/(4*geometry.radius_tree);
SWRabsDir_t.SWRabsTree	=	SWRabsDir_t.SWRabsTree*4*geometry.radius_tree*pi/(4*geometry.radius_tree);
SWRabsDiff_t.SWRabsTree	=	SWRabsDiff_t.SWRabsTree*4*geometry.radius_tree*pi/(4*geometry.radius_tree);

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
[G1WallSun,G2WallSun,dsWallSun]=conductive_heat_functions.Impervious_Conductive_Heat...
	(TemperatureC,TempVec,Anthropogenic,ParThermalWall,WallLayers,ParCalculation,itt,type);

% Conductive heat flux shaded wall
type	=	0;
[G1WallShade,G2WallShade,dsWallShade]=conductive_heat_functions.Impervious_Conductive_Heat...
	(TemperatureC,TempVec,Anthropogenic,ParThermalWall,WallLayers,ParCalculation,itt,type);

% Conductive heat flux impervious ground
[G1GroundImp,Tdp_ground_imp]=conductive_heat_functions.ForceRestore_conductive_heat_imp(TemperatureC,TempDamp,TempVec,...
	ParCalculation,ParThermalGround,FractionsGround,itt);

% Conductive heat flux bare ground
type	=	0; % bare
[G1GroundBare,Tdp_ground_bare]=conductive_heat_functions.ForceRestore_conductive_heat_soil(TemperatureC,TempDamp,Owater,TempVec,...
	ParCalculation,ParSoilGround,ParVegGround,ParVegTree,FractionsGround,itt,type);

% Conductive heat flux vegetated ground
type	=	1; % vegetated
[G1GroundVeg,Tdp_ground_veg]=conductive_heat_functions.ForceRestore_conductive_heat_soil(TemperatureC,TempDamp,Owater,TempVec,...
	ParCalculation,ParSoilGround,ParVegGround,ParVegTree,FractionsGround,itt,type);


%% Sensible and latent heat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[HfluxCanyon,LEfluxCanyon,ra_canyon,HumidityCan]=turbulent_heat_function.HeatFlux_canyon...
	(TemperatureC,Gemeotry_m,MeteoData,ParVegTree,ParTree);

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
	=turbulent_heat_function.HeatFlux_ground(TemperatureC,MeteoData,Gemeotry_m,geometry,FractionsGround,ParTree,...
	ParVegGround,ParVegTree,ParSoilGround,SoilPotW,Owater,Vwater,ExWater,Int,CiCO2Leaf,...
	ParInterceptionTree,ParCalculation,itt,SWRabsDir_t.SWRabsTree,SWRabsDiff_t.SWRabsTree,SWRabsDir_t.SWRabsGroundVeg,SWRabsDiff_t.SWRabsGroundVeg);

% Turbulent heat fluxes from sunlit and shaded wall to canyon
[HfluxWallSun,HfluxWallShade,Ewsun,Ewshade,LEwsun,LEwshade,RES_w1,RES_w2,rap_Zp1_In,rap_Zp2_In,...
	Hwsun1,Hwshade1,Hwsun2,Hwshade2,cp_atm,rho_atm,L_heat,Zp1,Zp2,rap_Zp1,rap_Zp2]...
	=turbulent_heat_function.HeatFlux_wall(TemperatureC,Gemeotry_m,MeteoData,ParVegTree,ParTree,ParVegGround,FractionsGround);


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
	Ycanyon(1)	=	TemperatureC(1,1)-303.15;
end
if FractionsGround.fbare>0
	Ycanyon(2)	=	SWRabs_t.SWRabsGroundBare + LWRabs_t.LWRabsGroundBare - G1GroundBare - HfluxGroundBare - LEfluxGroundBarePond - LEfluxGroundBareSoil;
else
	Ycanyon(2)	=	TemperatureC(1,2)-303.15;
end
if FractionsGround.fveg>0
	Ycanyon(3)	=	SWRabs_t.SWRabsGroundVeg + LWRabs_t.LWRabsGroundVeg - G1GroundVeg - HfluxGroundVeg - LEfluxGroundVegInt - LEfluxGroundVegPond - LEfluxGroundVegSoil - LTEfluxGroundVeg;
else
	Ycanyon(3)	=	TemperatureC(1,3)-303.15;
end
Ycanyon(4)	=	SWRabs_t.SWRabsWallSun + LWRabs_t.LWRabsWallSun - G1WallSun - HfluxWallSun;
Ycanyon(5)	=	SWRabs_t.SWRabsWallShade + LWRabs_t.LWRabsWallShade - G1WallShade - HfluxWallShade;
if ParTree.trees>0
	Ycanyon(6)	=	SWRabs_t.SWRabsTree + LWRabs_t.LWRabsTree - HfluxTree- LEfluxTreeInt - LTEfluxTree;		
else
	Ycanyon(6)	=	TemperatureC(1,6)-303.15;
end
Ycanyon(7)	=	G1WallSun - G2WallSun - dsWallSun; % Energy budget interior of sunlit wall
Ycanyon(8)	=	G1WallShade - G2WallShade - dsWallShade; % Energy budget interior of shaded wall
Ycanyon(9)	=	Anthropogenic.Qf_canyon + Cimp*FractionsGround.fimp*HfluxGroundImp + Cbare*FractionsGround.fbare*HfluxGroundBare + Cveg*FractionsGround.fveg*HfluxGroundVeg + geometry.hcanyon*HfluxWallSun + geometry.hcanyon*HfluxWallShade + Ctree*4*geometry.radius_tree*HfluxTree - HfluxCanyon;
Ycanyon(10)	=	Cimp*FractionsGround.fimp*LEfluxGroundImp + Cbare*FractionsGround.fbare*(LEfluxGroundBarePond+LEfluxGroundBareSoil) + ...
				Cveg*FractionsGround.fveg*(LEfluxGroundVegInt + LEfluxGroundVegPond + LEfluxGroundVegSoil + LTEfluxGroundVeg) + Ctree*4*geometry.radius_tree*(LEfluxTreeInt + LTEfluxTree) - LEfluxCanyon;

			
			
