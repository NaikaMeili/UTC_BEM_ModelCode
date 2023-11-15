function[dS]=HeatStorageChangeInternalMass(Tintmass,Tintmasstm1,ParThermalBuildFloor,Gemeotry_m,ParCalculation)

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

Ts			    =	Tintmass;
Ts_tm1	        =	Tintmasstm1;
dts			    =	ParCalculation.dts;
cv_floor        =	ParThermalBuildFloor.cv_floor_IntMass;
cv_wall         =	ParThermalBuildFloor.cv_wall_IntMass;
dzFloor         =   ParThermalBuildFloor.dzFloor;
dzWall          =   ParThermalBuildFloor.dzWall;
FloorHeight     =   ParThermalBuildFloor.FloorHeight;
BuildingHeight  =   Gemeotry_m.Height_canyon;
BuildingWidth   =   Gemeotry_m.Width_roof;

% Calculate number of floor layers to account for in internal mass
NrOfFloorLayers = max(round(BuildingHeight./FloorHeight)-1,0);

% Calculate internal mass for heat storage based on average floor thickness
% and number of floors and walls
IntFloorVolume      =   NrOfFloorLayers*dzFloor*BuildingWidth;
IntWallVolume       =   dzWall*BuildingHeight;

% This has to be normalized by building height as the shortwave, longwave
% and sensible heat exchange is calculate per m^2 of wall.
IntTotal_dz     = (IntFloorVolume+IntWallVolume)./BuildingHeight; 

% Calculation of average heat capacity of internal mass according to
% volumne of walls and ground
cv_int_total    = cv_floor.*IntFloorVolume/(IntFloorVolume+IntWallVolume) + cv_wall.*IntWallVolume/(IntFloorVolume+IntWallVolume);

%%%%%%%%%%%%%% COMPUTATION
dS  =   cv_int_total*IntTotal_dz/dts*(Ts-Ts_tm1);




