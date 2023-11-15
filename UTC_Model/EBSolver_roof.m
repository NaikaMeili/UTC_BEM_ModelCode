function[Yroof,G2Roof]=EBSolver_roof(TemperatureR,TemperatureB,TempVec_ittm,MeteoData,...
				Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,CiCO2Leaf_ittm,...
				Gemeotry_m,FractionsRoof,ParSoilRoof,PropOpticalRoof,ParThermalRoof,ParVegRoof,...
				HumidityAtm,Anthropogenic,ParCalculation,BEM_on,RESPreCalc,rsRoofPreCalc)


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
	=turbulent_heat_function.HeatFlux_roof(TemperatureR,TempVec_ittm,MeteoData,HumidityAtm,ParVegRoof,FractionsRoof,Gemeotry_m,...
	ParSoilRoof,ParCalculation,SoilPotW_ittm,Owater_ittm,Vwater_ittm,ExWater_ittm,Int_ittm,CiCO2Leaf_ittm,...
	SWRabs_dir_veg,SWRabs_diff_veg,RESPreCalc,rsRoofPreCalc);


%% Conductive heat fluxes rof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impervious conductive heat flux
[G1_roofimp,G2_roofimp,dS_roofimp]=conductive_heat_functions.ConductiveHeatFlux_RoofImp...
	(TemperatureR,TemperatureB,TempVec_ittm,Anthropogenic,ParThermalRoof,ParSoilRoof,ParCalculation,BEM_on);

% Conductive heat flux of green roof
[G1_roofveg,G2_roofveg,dS_roofveg]=conductive_heat_functions.ConductiveHeatFlux_GreenRoof...
	(TemperatureR,TemperatureB,TempVec_ittm,Anthropogenic,Owater_ittm,ParVegRoof,ParSoilRoof,ParThermalRoof,ParCalculation,BEM_on);

G2Roof = FractionsRoof.fimp.*G2_roofimp + FractionsRoof.fveg.*G2_roofveg;


%% Energy balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FractionsRoof.fimp>0
	Yroof(1)	=	SWR_abs_roofimp+LWR_abs_roofimp-Hroof_imp-G1_roofimp-LEroof_imp;  % Energy budget impervious roof
    Yroof(3)	=	G1_roofimp-G2_roofimp-dS_roofimp; % Energy budget concrete mass roof
else
	Yroof(1)	=	TemperatureR(1)-273.15;
    Yroof(3)	=	TemperatureR(3)-273.15;
end

if FractionsRoof.fveg>0
	Yroof(2)	=	SWR_abs_roofveg+LWR_abs_roofveg-Hroof_veg-G1_roofveg-LEroof_veg-LEroof_ground-LEroof_soil-LTEroof_veg;  % Energy budget vegetated roof
    Yroof(4)	=	G1_roofveg-G2_roofveg-dS_roofveg; % Energy budget concrete mass roof
else
	Yroof(2)	=	TemperatureR(2)-273.15;
    Yroof(4)	=	TemperatureR(4)-273.15;
end



end
