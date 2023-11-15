function[G1,G2,dS]=ConductiveHeatFlux_RoofImp(TemperatureR,TemperatureB,TempVec_ittm,Anthropogenic,ParThermalRoof,ParSoilRoof,ParCalculation,BEM_on)

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

Ts			=	TemperatureR(1,1);
Tint		=	TemperatureR(1,3);
if BEM_on ==1
    Tb      =	TemperatureB(1); % Ceiling temperature
else
    Tb      =   Anthropogenic.Tb;
end
Tint_tm1	=	TempVec_ittm.TRoofIntImp;
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

