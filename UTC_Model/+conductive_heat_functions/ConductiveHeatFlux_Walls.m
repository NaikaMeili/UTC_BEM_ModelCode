function[G1,G2,dS]=ConductiveHeatFlux_Walls(TemperatureC,TemperatureB,TempVec_ittm,TempVecB_ittm,Anthropogenic,...
                                            ParThermalWall,WallLayers,ParCalculation,type,ParWindows,BEM_on)

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
if BEM_on ==1
    Tb  =	TemperatureB(2);
else
    Tb  =	Anthropogenic.Tb;
end
Tint		=	TemperatureC(1,7);
Ts_tm1      =   TempVec_ittm.TWallSun;      
Tb_tm1      =   TempVecB_ittm.Tinwallsun;
Tint_tm1	=	TempVec_ittm.TWallIntSun;
lan_dry1	=	ParThermalWall.lan_dry;
lan_dry2	=	ParThermalWall.lan_dry;
dz1			=	WallLayers.dz1_wall;
dz2			=	WallLayers.dz2_wall;
cv_s1		=	ParThermalWall.cv_s;
cv_s2		=	ParThermalWall.cv_s;
dts			=	ParCalculation.dts;
elseif type == 0 % shade
Ts			=	TemperatureC(1,5);

if BEM_on==1
    Tb  =	TemperatureB(3);
else
    Tb  =	Anthropogenic.Tb;
end
Tint		=	TemperatureC(1,8);
Ts_tm1      =   TempVec_ittm.TWallShade;      
Tb_tm1      =   TempVecB_ittm.Tinwallshd;
Tint_tm1	=	TempVec_ittm.TWallIntShade;
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


UvalueWindow    = ParWindows.Uvalue;  % W/m^2K, U-value = Thermal conductivity / thickness
cv_glass        = ParWindows.cv_glass;
dzglass         = ParWindows.dztot;
RatioWindowToWall = ParWindows.GlazingRatio;
if RatioWindowToWall>0.95
    RatioWindowToWall = 0.95;
end


%%%%%%%%%%%%%% COMPUTATION
% Conductive heat flux through walls
%---------------------------------------------------
G1wall      =	lan_dry1*(Ts-Tint)/dz1;
G2wall      =	lan_dry2*(Tint-Tb)/dz2;
dSwall      =	(cv_s1+cv_s2)/2*(dz1+dz2)/dts*(Tint-Tint_tm1);

% Conductive heat flux through windows, assume absorption is negligible
%---------------------------------------------------
G1window    =   UvalueWindow*(Ts-Tb);
G2window    =   UvalueWindow*(Ts-Tb);
dSwindow    =   cv_glass*dzglass/dts*((Ts+Tb)/2 - (Ts_tm1+Tb_tm1)/2);

% Average conductive heat fluxes according to window to wall ratio
if BEM_on==1
    G1 = (1-RatioWindowToWall)*G1wall + RatioWindowToWall*G1window;
    G2 = (1-RatioWindowToWall)*G2wall + RatioWindowToWall*G2window;
    dS = (1-RatioWindowToWall)*dSwall + RatioWindowToWall*dSwindow;
else
    G1 = G1wall;
    G2 = G2wall;
    dS = dSwall;
end


