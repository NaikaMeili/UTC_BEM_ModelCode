function[G1,G2,dS]=Impervious_Conductive_HeatRoof(TemperatureR,TempVec,Anthropogenic,ParThermalRoof,ParSoilRoof,ParCalculation,itt)

%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dts		=	time step [s] = 1 h = 3600 s
% dz		=	depth of layer [m]
% Ts		=	surface temperature [K]
% Tb		=	Constant building interior temperature [K]
% Tint		=	Temperature of concrete [K]
% lan_dry	=	Thermal conductivity dry solid [W/m K]
% cv_s		=	Volumetric heat capacity solid [J/m^3 K]

%%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% G1		=	Heat flux from surface to concrete interior [W/m^2]
% G2		=	Heat flux from concrete interior to building interior [W/m^2]
% dS		=	Energy storage in the roof

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

Ts			=	TemperatureR(1,1);
Tint		=	TemperatureR(1,3);
Tb			=	Anthropogenic.Tb;
Tint_tm1	=	TempVec.TRoofIntImp(itt,1);
lan_dry1	=	ParThermalRoof.lan_dry_imp;
lan_dry2	=	ParThermalRoof.lan_dry_imp;
dz1			=	ParSoilRoof.dz1;
dz2			=	ParSoilRoof.dz2;
cv_s1		=	ParThermalRoof.cv_s_imp;
cv_s2		=	ParThermalRoof.cv_s_imp;
dts			=	ParCalculation.dts;

	

%%%%%%%%%%%%%% COMPUTATION
G1			=	lan_dry1*(Ts-Tint)/dz1; % Soil Heat Flux [W/m^2];
G2			=	lan_dry2*(Tint-Tb)/dz2;
dS			=	(cv_s1+cv_s2)/2*(dz1+dz2)/dts*(Tint-Tint_tm1);

end
