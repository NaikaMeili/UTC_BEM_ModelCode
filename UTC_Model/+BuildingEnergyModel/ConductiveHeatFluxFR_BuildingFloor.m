function[G,Tdp]=ConductiveHeatFluxFR_BuildingFloor(Tinground,TingroundDamptm1,Tingroundtm1,...
	ParCalculation,ParThermalBulidFloor)


Ts				=	Tinground;	
Tdptm1			=	TingroundDamptm1;
Tstm1			=	Tingroundtm1;
dts				=	ParCalculation.dts;
lan_dry_imp		=	ParThermalBulidFloor.lan_ground_floor;
cv_s_imp		=	ParThermalBulidFloor.cv_ground_floor;


tau		=	86400; %% [s] time constant
CTt		=	2*(sqrt(pi/(lan_dry_imp*cv_s_imp*tau))); %%  [K m^2/J] Total Thermal Capacity Soil

[G,Tdp]	=	conductive_heat_functions.Soil_Heat(dts,Ts-273.15,Tstm1-273.15,Tdptm1-273.15,CTt);
Tdp		=	Tdp + 273.15;
