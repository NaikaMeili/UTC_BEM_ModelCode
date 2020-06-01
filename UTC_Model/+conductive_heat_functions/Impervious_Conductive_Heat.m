function[G1,G2,dS]=Impervious_Conductive_Heat(TemperatureC,TempVec,Anthropogenic,ParThermalWall,WallLayers,ParCalculation,itt,type)

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

if type == 1 % sun
Ts			=	TemperatureC(1,4);
Tb			=	Anthropogenic.Tb;
Tint		=	TemperatureC(1,7);
Tint_tm1	=	TempVec.TWallIntSun(itt,1);
lan_dry1	=	ParThermalWall.lan_dry;
lan_dry2	=	ParThermalWall.lan_dry;
dz1			=	WallLayers.dz1_wall;
dz2			=	WallLayers.dz2_wall;
cv_s1		=	ParThermalWall.cv_s;
cv_s2		=	ParThermalWall.cv_s;
dts			=	ParCalculation.dts;
elseif type == 0 % shade
Ts			=	TemperatureC(1,5);
Tb			=	Anthropogenic.Tb;
Tint		=	TemperatureC(1,8);
Tint_tm1	=	TempVec.TWallIntShade(itt,1);
lan_dry1	=	ParThermalWall.lan_dry;
lan_dry2	=	ParThermalWall.lan_dry;
dz1			=	WallLayers.dz1_wall;
dz2			=	WallLayers.dz2_wall;
cv_s1		=	ParThermalWall.cv_s;
cv_s2		=	ParThermalWall.cv_s;
dts			=	ParCalculation.dts;
else
	disp('please, enter sun or shade for sunlit or shaded wall')
end
	

%%%%%%%%%%%%%% COMPUTATION
G1			=	lan_dry1*(Ts-Tint)/dz1; % Soil Heat Flux [W/m^2];
G2			=	lan_dry2*(Tint-Tb)/dz2;
dS			=	(cv_s1+cv_s2)/2*(dz1+dz2)/dts*(Tint-Tint_tm1);

end
