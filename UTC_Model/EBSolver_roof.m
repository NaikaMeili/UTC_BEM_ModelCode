function[Yroof]=EBSolver_roof(TemperatureR,TempVec,MeteoData,itt_tm1,...
				Int,ExWater,Vwater,Owater,SoilPotW,CiCO2Leaf,...
				Gemeotry_m,FractionsRoof,ParSoilRoof,PropOpticalRoof,ParThermalRoof,ParVegRoof,...
				HumidityAtm,Anthropogenic,ParCalculation)

% TemperatureR(1,1) = Troof_imp
% TemperatureR(1,2) = Troof_veg
% TemperatureR(1,3) = Troof_interior_imp
% TemperatureR(1,4) = Troof_interior_veg


%% Shortwave radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SWRabs_dir_veg		=	(1-PropOpticalRoof.aveg)*MeteoData.SW_dir;	% Absorbed direct shortwave radiation by vegetated roof [W/m^2]
SWRabs_diff_veg		=	(1-PropOpticalRoof.aveg)*MeteoData.SW_diff;	% Absorbed diffuse shortwave radiation by vegetated roof [W/m^2]
SWR_abs_roofveg		=	(1-PropOpticalRoof.aveg)*(MeteoData.SW_dir+MeteoData.SW_diff);	% Absorbed total shortwave radiation by vegetated roof [W/m^2]
SWR_abs_roofimp		=	(1-PropOpticalRoof.aimp)*(MeteoData.SW_dir+MeteoData.SW_diff);	% Absorbed total shortwave radiation by impervious roof [W/m^2]

SWR_out_roofveg		=	PropOpticalRoof.aveg*(MeteoData.SW_dir+MeteoData.SW_diff);
SWR_out_roofimp		=	PropOpticalRoof.aimp*(MeteoData.SW_dir+MeteoData.SW_diff);

%% Longwave radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bolzm				=	5.67*10^(-8);	% Stefan-Boltzmann constant [W*m^-2*K-4], Lecture notes hydrology 2
LWR_abs_roofveg		=	MeteoData.LWR-(PropOpticalRoof.eveg*bolzm*(TemperatureR(1,2))^4+(1-PropOpticalRoof.eveg)*MeteoData.LWR);	% Total absorbed longwave radiation by vegetated roof [W/m^2]
LWR_abs_roofimp		=	MeteoData.LWR-(PropOpticalRoof.eimp*bolzm*(TemperatureR(1,1))^4+(1-PropOpticalRoof.eimp)*MeteoData.LWR);	% Total absorbed longwave radiation by impervious roof [W/m^2]

LWR_out_roofveg		=	PropOpticalRoof.eveg*bolzm*(TemperatureR(1,2))^4+(1-PropOpticalRoof.eveg)*MeteoData.LWR;
LWR_out_roofimp		=	PropOpticalRoof.eimp*bolzm*(TemperatureR(1,1))^4+(1-PropOpticalRoof.eimp)*MeteoData.LWR;

%% Sensible and latent heat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Hroof_imp,Hroof_veg,Eroof_imp,Eroof_veg,Eroof_ground,Eroof_soil,TEroof_veg,...
	LEroof_imp,LEroof_veg,LEroof_ground,LEroof_soil,LTEroof_veg,...
	Ci_sun_roof,Ci_shd_roof,ra,rb_L,rap_L,r_soil,rs_sun,rs_shd]...
	=turbulent_heat_function.HeatFlux_roof(TemperatureR,MeteoData,HumidityAtm,ParVegRoof,FractionsRoof,Gemeotry_m,...
	ParSoilRoof,ParCalculation,SoilPotW,Owater,Vwater,ExWater,Int,CiCO2Leaf,...
	itt_tm1,SWRabs_dir_veg,SWRabs_diff_veg);


%% Ground heat flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impervious conductive heat flux
[G1_roofimp,G2_roofimp,dS_roofimp]=conductive_heat_functions.Impervious_Conductive_HeatRoof...
	(TemperatureR,TempVec,Anthropogenic,ParThermalRoof,ParSoilRoof,ParCalculation,itt_tm1);

% Vegetated ground heat flux
[G1_roofveg,G2_roofveg,dS_roofveg]=conductive_heat_functions.Soil_Conductive_Heat...
	(TemperatureR,TempVec,Anthropogenic,Owater,ParVegRoof,ParSoilRoof,ParThermalRoof,ParCalculation,itt_tm1);


%% Energy balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Yroof(1)	=	SWR_abs_roofimp+LWR_abs_roofimp-Hroof_imp-G1_roofimp-LEroof_imp;  % Energy budget impervious roof
Yroof(2)	=	SWR_abs_roofveg+LWR_abs_roofveg-Hroof_veg-G1_roofveg-LEroof_veg-LEroof_ground-LEroof_soil-LTEroof_veg;  % Energy budget vegetated roof
Yroof(3)	=	G1_roofimp-G2_roofimp-dS_roofimp; % Energy budget concrete mass roof
Yroof(4)	=	G1_roofveg-G2_roofveg-dS_roofveg; % Energy budget concrete mass roof


end
