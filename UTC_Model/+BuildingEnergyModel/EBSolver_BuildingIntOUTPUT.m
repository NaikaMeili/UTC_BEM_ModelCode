function[HbuildInt,LEbuildInt,GbuildInt,SWRabsB,LWRabsB,Tdpfloor,WasteHeat,EnergyUse]=EBSolver_BuildingIntOUTPUT(TemperatureB,MeteoData,SWRinWsun,SWRinWshd,...
    TempVecB_ittm,Tairout,qairout,G2Roof,G2WallSun,G2WallShade,TBdamp_ittm,...
    Gemeotry_m,PropOpticalIndoors,ParHVAC,ParCalculation,ParThermalBulidFloor,...
    TempVec,G1Roof,G1WallSun,G1WallShade,dsWallSun,dsWallShade)

% Simple Building energy model
%--------------------------------------------------------------------------

% Temperature vector:
% TemperatureB(:,1)		=	Temperature ceiling
% TemperatureB(:,2)		=	Temperature sunlit wall
% TemperatureB(:,3)		=	Temperature shaded wall
% TemperatureB(:,4)		=	Temperature ground
% TemperatureB(:,5)		=	Temperature air
% TemperatureB(:,6)		=	Humidity air


% Internal temperatures and humidity that are calculated
Tceiling    = TemperatureB(1);
Tinwallsun  = TemperatureB(2);
Tinwallshd  = TemperatureB(3);
Tinground   = TemperatureB(4);
Tbin        = TemperatureB(5);
qbin        = TemperatureB(6);

% Assign variable
Tbin_tm1 = TempVecB_ittm.Tbin;
qbin_tm1 = TempVecB_ittm.qbin;
Tingroundtm1 = TempVecB_ittm.Tinground; 
TingroundDamptm1    =   TBdamp_ittm.TingroundDamp;
Tatm     = MeteoData.Tatm;
Pre      = MeteoData.Pre;
ea       = MeteoData.ea;


% Urban geometry (in metres)
Hbuild = Gemeotry_m.Height_canyon;
Wroof = Gemeotry_m.Width_roof;

% albedo/emissivities
abc = PropOpticalIndoors.abc; abw = PropOpticalIndoors.abw; abg = PropOpticalIndoors.abg;
ec = PropOpticalIndoors.ec; eg = PropOpticalIndoors.eg; ew = PropOpticalIndoors.ew;

% Aircon parameters
ACon = ParHVAC.ACon;
ACH = ParHVAC.ACH; % air changes per hour (1/h)
Tsetpoint = ParHVAC.Tsetpoint;
COP = ParHVAC.COP;

% Temporal resolution
dth = ParCalculation.dth; % (h)
dts = ParCalculation.dts; % (s) number of seconds in 1 hour


if ACon==1
    Tbin_tm1 = Tsetpoint;
    Tbin = Tsetpoint;
    TemperatureB(5) = Tsetpoint;
end


% Constants, calculated as function of air temperature, pressure, humidity 
%--------------------------------------------------------------------------
L_heat  =	1000*(2501.3 - 2.361*(Tatm-273.15)); % Latent heat vaporization/condensation [J/kg], hwv=1000.*2501;% J/kg Specific enthalpy of water vapour
Cpa     =	1005+(((Tatm-273.15)+23.15)^2)/3364; % Specific heat capacity of the air [J/kg K]
Cpwv	=   1000.*1.84;   % J/kg degC Specific heat capacity of water vapour
rho_atm =	(Pre/(287.04*Tatm))*(1-(ea/Pre)*(1-0.622));	% dry air density at atmosphere [kg/m^3]


% Interior building surface fluxes
%--------------------------------------------------------------------------
% Shortwave radiation 
[SWRinB,SWRoutB,SWRabsB,SWREBB]=BuildingEnergyModel.SWRabsIndoors(SWRinWsun,SWRinWshd,Hbuild,Wroof,abc,abw,abg);

% Longwave radiation
% Windows are assumed opaque for longwave radiation
[LWRinB,LWRoutB,LWRabsB,LWREBB]=BuildingEnergyModel.LWRabsIndoors(Tinwallsun,Tinwallshd,Tceiling,Tinground,Hbuild,Wroof,ec,eg,ew);

% Sensible heat fluxes, positive flux indicates flux from surface to air
[HbinWallSun,HbinWallshd,HBinRoof,HBinGround]=...
    BuildingEnergyModel.SensibleHeatFluxBuildingInterior(Tbin,Tinwallsun,Tinwallshd,Tceiling,Tinground);

HbuildIn =  Hbuild.*(HbinWallSun+HbinWallshd) + Wroof.*HBinRoof + ...
            Wroof.*HBinGround; % (W/m) Building width

% Conductive heat fluxes
[Gfloor,Tdpfloor]=BuildingEnergyModel.ForceRestore_conductive_heat_BuildFloor(Tinground,TingroundDamptm1,Tingroundtm1,...
	ParCalculation,ParThermalBulidFloor);


%--------------------------------------------------------------------------
% Internal sensible and latent heat gain (W/m), positive flux indicates
% added to indoor air
[Hequip,Hpeople,LEequip,LEpeople]=BuildingEnergyModel.IndoorSensibleLatentHeatSource();


% Sensible heat load through ventilation, 
% Positive flux indicates that outdoor air is warmer than indoor air, 
% i.e. heat is added from outdoors to indoors
% Negative flux indicates that air is colder outdoors than indoors and
% hence, heat is removed from the indoor air.
Vbuild = Wroof.*Hbuild;
Hvent   =   (ACH.*dth.*Vbuild)./3600.*Cpa.*rho_atm.*(Tairout-Tbin); % (1/h)*(h)*m^2*J/(kg K)*kg/m^3*K, % (W/m) Building length
LEvent  =   (ACH.*dth.*Vbuild)./3600.*rho_atm.*L_heat*(qairout-qbin); % (1/h)*(h)*m^2*(kg/m^3)*(J/kg)*(kg/kg), % (W/m) Building length

% AC modlue
%--------------------------------
if ACon==1
    H_AC    = HbuildIn + Hvent + Hequip + Hpeople;
    LE_AC   = LEvent + LEequip + LEpeople;
else
    H_AC    = 0;
    LE_AC   = 0;
end

H_air = Vbuild.*Cpa.*rho_atm.*(Tbin-Tbin_tm1)./dts; % (W/m), m^2*J/(kg K)*(kg/m^3)*K/s
LE_air = Vbuild.*rho_atm.*L_heat*(qbin-qbin_tm1)./dts; % (W/m), m^2*(kg/m^3)*(J/kg)*(kg/kg)/s


% Energy balance for individual surfaces
YBuildInt(1)	=	SWRabsB.SWRabsCeiling + LWRabsB.LWRabsCeiling + G2Roof - HBinRoof;
YBuildInt(2)	=	SWRabsB.SWRabsWallsun + LWRabsB.LWRabsWallsun + G2WallSun - HbinWallSun;
YBuildInt(3)	=	SWRabsB.SWRabsWallshade + LWRabsB.LWRabsWallshade + G2WallShade - HbinWallshd;
YBuildInt(4)	=	SWRabsB.SWRabsGround + LWRabsB.LWRabsGround - Gfloor - HBinGround;
YBuildInt(5)    =   HbuildIn + Hvent + Hequip + Hpeople - H_AC - H_air; % (W/m)
YBuildInt(6)    =   LEvent + LEequip + LEpeople - LE_AC - LE_air; % (W/m)

% [G1Roof,G1WallSun,G1WallShade,0;...
%     G2Roof,G2WallSun,G2WallShade,Gfloor]
% 
% [dsWallSun,dsWallShade]
% 
% [TempVec(4),TempVec(7),TempVec(5),TempVec(8)]-273.15
% TempVec-273.15
% 
% TemperatureB-273.15


% Building energy consumption
EnergyUse   = (H_AC+LE_AC)./COP; % W, COP for heating and cooling could be different
WasteHeat   = (COP+1)./COP.*(H_AC+LE_AC); % W


%--------------------------------------------------------------------------
% Prepare energy flux outputs
% Sensible heat fluxes-----------------------------------------------------
HbuildInt.HBinRoof = HBinRoof;
HbuildInt.HbinWallSun = HbinWallSun;
HbuildInt.HbinWallshd = HbinWallshd;
HbuildInt.HBinGround = HBinGround;
HbuildInt.HbuildInSurf = HbuildIn;
HbuildInt.Hvent = Hvent;
HbuildInt.Hequip = Hequip;
HbuildInt.Hpeople = Hpeople;
HbuildInt.H_AC = H_AC;
HbuildInt.H_air = H_air;

% Latent heat fluxes-----------------------------------------------------
LEbuildInt.LEvent = LEvent;
LEbuildInt.LEequip = LEequip;
LEbuildInt.LEpeople = LEpeople;
LEbuildInt.LE_AC = LE_AC;
LEbuildInt.LE_air = LE_air;

% Conductive heat fluxes-----------------------------------------------------
GbuildInt.G2Roof = G2Roof;
GbuildInt.G2WallSun = G2WallSun;
GbuildInt.G2WallShade = G2WallShade;
GbuildInt.Gfloor = Gfloor;












