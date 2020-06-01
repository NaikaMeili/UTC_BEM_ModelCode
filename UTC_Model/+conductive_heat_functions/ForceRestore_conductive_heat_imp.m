function[G,Tdp]=ForceRestore_conductive_heat_imp(TemperatureC,TempDamp,TempVec,...
	ParCalculation,ParThermalGround,FractionsGround,itt)

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
Tdptm1			=	TempDamp.TDampGroundImp(itt,1);
Tstm1			=	TempVec.TGroundImp(itt,1);
dts				=	ParCalculation.dts;
lan_dry_imp		=	ParThermalGround.lan_dry_imp;
cv_s_imp		=	ParThermalGround.cv_s_imp;
Cimp			=	(FractionsGround.fimp>0);

tau		=	86400; %% [s] time constant
CTt		=	2*(sqrt(pi/(lan_dry_imp*cv_s_imp*tau))); %%  [K m^2/J] Total Thermal Capacity Soil

[G,Tdp]	=	conductive_heat_functions.Soil_Heat(dts,Ts-273.15,Tstm1-273.15,Tdptm1-273.15,CTt);
Tdp		=	Tdp + 273.15;
G		=	G*Cimp;
