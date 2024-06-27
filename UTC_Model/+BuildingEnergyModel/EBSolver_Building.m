function[YBuildInt,WasteHeat]=EBSolver_Building(TemperatureC,TemperatureB,TempVecB_ittm,TempVec_ittm,Humidity_ittm,MeteoData,...
    SWRinWsun,SWRinWshd,G2Roof,G2WallSun,G2WallShade,TempDamp_ittm,SWRabs_t,...
    Gemeotry_m,PropOpticalIndoors,ParHVAC,ParCalculation,ParThermalBulidingInt,ParWindows,BEM_on,HVACSchedule)


% Simple Building energy model
%--------------------------------------------------------------------------
% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wasteheat = sensible and latent heat emissions from AC and from venitlation into the canyon

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TemperatureC: outdoor temperatures of buildings and air and humidity [K]
% TemperatureB: Building internal temperature and air and humidity [K]
% TempVecB_ittm: Building internal temperature and humidity from previous time step
% TempDamp_ittm: Ground dampening temperature from previous time step
% MeteoData: Atmospheric forcing conditions
% Gemeotry_m: Urban geometry (in metres)
% PropOpticalIndoors: Building internal albedo/emissivities
% ParHVAC.ACon: Parameter specifying if AC is on (=1) or off (=0)
% ParHVAC.ACH: Ventilation rate: air changes per hour (1/h)
% ParHVAC.Tsetpoint: Setpoint temperature for air conditioning
% ParHVAC.COP: Coefficient of performance of air conditioning
% ParCalculation: Temporal resolution of calculation


% Initialization of prognostic variables in building interior
%--------------------------------------------------------------------------
% Internal temperatures and humidity that are calculated [K]
Tceiling    = TemperatureB(1); % Temperature ceiling
Tinwallsun  = TemperatureB(2); % Temperature sunlit wall
Tinwallshd  = TemperatureB(3); % Temperature shaded wall
Twindow     = TemperatureB(4); % Temperature windows
Tinground   = TemperatureB(5); % Temperature ground
Tintmass    = TemperatureB(6); % Temperature internal mass
Tbin        = TemperatureB(7); % Air temperature building interior
qbin        = TemperatureB(8); % Specific humidity air building interior [kg/kg]

% External canyon air temperature and specific humidity
Tairout = TemperatureC(9);  % Air temperaure at outdoor calculation height (hdisp + zom) [K]
qairout = TemperatureC(10); % Specific humidity at outdoor calucation height (hdisp + zom) [kg/kg]

% Urban geometry [m]
Hwall   = Gemeotry_m.Height_canyon;     % Building height [m]
Wroof   = Gemeotry_m.Width_roof;        % Roof width [m]
Wcan    = Gemeotry_m.Width_canyon;      % Street/canyon width [m]

% Temperature [K] and specific humidity [kg/kg] conditions from previous time step
Tbin_tm1            = TempVecB_ittm.Tbin;               % Building internal temperature from previous timestep
qbin_tm1            = TempVecB_ittm.qbin;               % Building internal humidity from previous timestep
Tingroundtm1        = TempVecB_ittm.Tinground;          % Ground temperature from previous timestep 
TingroundDamptm1    = TempDamp_ittm.TDampGroundBuild;   % Within ground dampening temperature from previous time step
Tintmasstm1         = TempVecB_ittm.Tintmass;           % Ground temperature from previous timestep   

Tsurftm1 = (Wroof*(TempVecB_ittm.Tinground + TempVecB_ittm.Tceiling) + Hwall*(TempVecB_ittm.Tinwallshd + TempVecB_ittm.Tinwallsun))/(2*Wroof + 2*Hwall);

% Atmospheric forcing conditions
Tatm     = MeteoData.Tatm;  % Air temperature at atmospheric forcing height Zatm [K]
Pre      = MeteoData.Pre;   % Air pressure at atmospheric forcing height Zatm [Pa]
ea       = MeteoData.ea;    % Vapor pressure at atmospheric forcing height Zatm [Pa]

% Albedo [-]
abc = PropOpticalIndoors.abc; % Ceiling
abw = PropOpticalIndoors.abw; % Internal wall 
abg = PropOpticalIndoors.abg; % Internal ground
abm = PropOpticalIndoors.abm; % Internal mass (i.e. internal walls)
% Emissivity [-]
ec = PropOpticalIndoors.ec; % Ceiling 
eg = PropOpticalIndoors.eg; % Internal wall 
ew = PropOpticalIndoors.ew; % Internal ground
em = PropOpticalIndoors.em; % Internal mass (i.e. internal walls)

% Heating and cooling parameters
AC_on           = ParHVAC.ACon;                     % Parameter specifying if AC is on (=1) or off (=0)
Heat_on         = ParHVAC.Heatingon;                % Parameter specifying if building internal heating is on and off (on = 1, off = 0)
TspCooling      = ParHVAC.TsetpointCooling;         % Cooling set-point temperature [K]
TspHeating      = ParHVAC.TsetpointHeating;         % Heating set-point temperature [K]
RHsptCooling    = ParHVAC.RHsetpointCooling./100;   % Cooling set-point humidity [-]
RHsptHeating    = NaN;                              % Heating set-point humidity, currently not used in the model
COP_AC          = ParHVAC.COPAC;                    % Coefficient of performance for AC [-]
COP_Heat        = ParHVAC.COPHeat;                  % Coefficient of performance for Heating [-]
ACH             = ParHVAC.ACH;                      % Ventilation rate: air changes per hour [1/h]
f_AC_LEToQ      = ParHVAC.f_ACLatentToQ;            % Fraction (0-1) of latent heat removed by AC that is condensed and ends up in the wastewater/ruonff.

AC_onCool       = ParHVAC.AC_onCool;                % Parameter specifying if AC is on (=1) or off (=0)
AC_onDehum      = ParHVAC.AC_onDehum;               % Parameter specifying if building internal heating is on and off (on = 1, off = 0)

% % Control in case heating and cooling set-points are the same
% if TspCooling == TspHeating  && AC_on==1 && Heat_on==1
%     TspHeating = TspHeating-0.01;
% elseif TspHeating > TspCooling && AC_on==1 && Heat_on==1
%     disp("Cooling set point is lower than heating setpoint, please double check")
% end

% Temporal resolution of calculation
dth = ParCalculation.dth; % (h)
dts = ParCalculation.dts; % (s)

% Calculate different humidity metrics for indoor humidity regulation
%--------------------------------------------------------------------------
esat_TspCooling  =    611*exp(17.27*(TspCooling-273.16)/(237.3+(TspCooling-273.16)));	% Saturation vapor pressure at set-point temperature [Pa]
ea_RHspCooling   =    RHsptCooling.*esat_TspCooling; % Indoor vapor pressure at building aircon set-point [Pa]
q_RHspCooling    =    0.622*ea_RHspCooling/(MeteoData.Pre-0.378*ea_RHspCooling);    % Specifc humidity of air at building aircon set-point [kg/kg]

% Constants, calculated as function of air temperature, pressure, humidity 
%--------------------------------------------------------------------------
L_heat  =	1000*(2501.3 - 2.361*(Tatm-273.15)); % Latent heat vaporization/condensation [J/kg], hwv=1000.*250 [J/kg] Specific enthalpy of water vapour
Cpa     =	1005+(((Tatm-273.15)+23.15)^2)/3364; % Specific heat capacity of the air [J/kg K]
Cpwv	=   1000.*1.84;   % Specific heat capacity of water vapour [J/kg K] 
rho_atm =	(Pre/(287.04*Tatm))*(1-(ea/Pre)*(1-0.622));	% Dry air density at atmospheric pressure [kg/m^3]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Building interior energy fluxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shortwave radiation [W/m^2] surface area of each element
%--------------------------------------------------------------------------
if ParThermalBulidingInt.IntMassOn == 1
    [SWRinB,SWRoutB,SWRabsB]=BuildingEnergyModel.SWRabsIndoors(SWRinWsun,SWRinWshd,Hwall,Wroof,abc,abw,abg,abm);
else
    [SWRinB,SWRoutB,SWRabsB,SWREBB]=BuildingEnergyModel.SWRabsIndoorsNoIntMass(SWRinWsun,SWRinWshd,Hwall,Wroof,abc,abw,abg);
    SWRinB.SWRinInternalMass = 0; SWRoutB.SWRoutInternalMass = 0; SWRabsB.SWRabsInternalMass = 0;
end

% Longwave radiation [W/m^2] surface area of each element
%--------------------------------------------------------------------------
% Windows are assumed opaque for longwave radiation
if ParThermalBulidingInt.IntMassOn == 1
    [LWRinB,LWRoutB,LWRabsB]=BuildingEnergyModel.LWRabsIndoors(Tinwallsun,Tinwallshd,Tceiling,Tinground,Tintmass,Hwall,Wroof,ec,eg,em,ew);
else
    [LWRinB,LWRoutB,LWRabsB,LWREBB]=BuildingEnergyModel.LWRabsIndoorsNoIntMass(Tinwallsun,Tinwallshd,Tceiling,Tinground,Hwall,Wroof,ec,eg,ew);
    LWRinB.LWRinInternalMass = 0; LWRoutB.LWRoutInternalMass = 0; LWRabsB.LWRabsInternalMass = 0;
end

% Sensible heat fluxes [W/m^2] surface area of each element 
% Positive flux indicates flux from surface to air
%--------------------------------------------------------------------------
[HbinWallSun,HbinWallshd,HBinRoof,HBinGround,HbinIntMass,HbinWindow]=...
    BuildingEnergyModel.SensibleHeatFluxBuildingInterior(Tbin,Tinwallsun,Tinwallshd,Tceiling,Tinground,Tintmass,Twindow);

% Calculate total internal sensible heat flux per m^2 building ground area
if ParThermalBulidingInt.IntMassOn == 1
    HbuildIn =  Hwall./Wroof.*(HbinWallSun+HbinWallshd) + HBinRoof + HBinGround + 2.*Hwall./Wroof.*HbinIntMass; % (W/m^2) Building ground area
else
    HbuildIn =  Hwall./Wroof.*(HbinWallSun+HbinWallshd) + HBinRoof + HBinGround;  % (W/m^2) Building ground area
    HbinIntMass = 0; 
end

% Conductive heat flux at building ground/floor [W/m^2] ground area
%--------------------------------------------------------------------------
[Gfloor,Tdpfloor]=BuildingEnergyModel.ConductiveHeatFluxFR_BuildingFloor(Tinground,TingroundDamptm1,Tingroundtm1,...
	ParCalculation,ParThermalBulidingInt);

% Heat storage of internal mass [W/m^2] wall height area
%--------------------------------------------------------------------------
if ParThermalBulidingInt.IntMassOn == 1
    [dSinternalMass]=BuildingEnergyModel.HeatStorageChangeInternalMass(Tintmass,Tintmasstm1,ParThermalBulidingInt,Gemeotry_m,ParCalculation);
else
    dSinternalMass = 0;
end

% Internal sensible and latent heat sources (W/m^2) per ground area, positive flux indicates added to indoor air
%--------------------------------------------------------------------------
%[Hequip,Hpeople,LEequip,LEpeople]=BuildingEnergyModel.IndoorSensibleLatentHeatSource();
Hequip      = HVACSchedule.Hequip;
Hpeople     = HVACSchedule.Hpeople;
LEequip     = HVACSchedule.LEequip;
LEpeople    = HVACSchedule.LEpeople;


% Sensible and latent heat load due to ventilation (W/m^2) per ground area, (air exchange between indoor and outdoor air) 
% Positive flux indicates that outdoor air is warmer than indoor air and heat is added from outdoors to indoors
% Negative flux indicates that air is colder outdoors than indoors and heat is removed from the indoor air.
%--------------------------------------------------------------------------
Vbuild  =   Wroof.*Hwall;
Hvent   =   (ACH.*dth.*Vbuild)./3600.*Cpa.*rho_atm.*(Tairout-Tbin)./Wroof; % (1/h)*(h)*m^2*J/(kg K)*kg/m^3*K, % (W/m2) Ground area
LEvent  =   (ACH.*dth.*Vbuild)./3600.*rho_atm.*L_heat*(qairout-qbin)./Wroof; % (1/h)*(h)*m^2*(kg/m^3)*(J/kg)*(kg/kg), % (W/m2) Ground area


% Change in heat storage and humidity storage in indoor air (W/m^2) per ground area
%--------------------------------------------------------------------------
dSH_air = Vbuild.*Cpa.*rho_atm.*(Tbin-Tbin_tm1)./dts./Wroof; % (W/m^2), m^2*J/(kg K)*(kg/m^3)*K/s
dSLE_air = Vbuild.*rho_atm.*L_heat*(qbin-qbin_tm1)./dts./Wroof; % (W/m^2), m^2*(kg/m^3)*(J/kg)*(kg/kg)/s


% AC and heating module
%--------------------------------------------------------------------------
[AC_on,AC_onCool,AC_onDehum,Heat_on,H_AC_Heat,LE_AC_Heat]=...
    BuildingEnergyModel.AC_HeatingModule(AC_on,Heat_on,AC_onCool,AC_onDehum,ParHVAC,...
    HbuildIn,Hvent,Hequip,Hpeople,dSH_air,LEvent,LEequip,LEpeople,dSLE_air);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy balance for individual surfaces and air volumne
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if BEM_on ==1
    % Ceiling energy balance
    YBuildInt(1)	=	SWRabsB.SWRabsCeiling + LWRabsB.LWRabsCeiling + G2Roof - HBinRoof;
    % Sunlit wall energy balance
    YBuildInt(2)	=	SWRabsB.SWRabsWallsun + LWRabsB.LWRabsWallsun + G2WallSun - HbinWallSun;
    % Shaded wall energy balance
    YBuildInt(3)	=	SWRabsB.SWRabsWallshd + LWRabsB.LWRabsWallshd + G2WallShade - HbinWallshd;
    % Window energy balance, currently not used in the model
    YBuildInt(4)	=	TemperatureB(4) - 273.15;
    % Ground/floor energy balance
    YBuildInt(5)	=	SWRabsB.SWRabsGround + LWRabsB.LWRabsGround - Gfloor - HBinGround;
    % Internal mass energy balance
    if ParThermalBulidingInt.IntMassOn == 1
    YBuildInt(6)    =   SWRabsB.SWRabsInternalMass + LWRabsB.LWRabsInternalMass - 2*HbinIntMass - dSinternalMass;
    else
    YBuildInt(6)    =   Tintmass - 273.15;
    end
    YBuildInt(7)    =   HbuildIn + Hvent + Hequip + Hpeople - dSH_air  - H_AC_Heat + ...
                        AC_onCool*1000*(Tbin - TspCooling) + Heat_on*1000*(Tbin - TspHeating); % (W/m2) building ground area
    YBuildInt(8)    =   LEvent + LEequip + LEpeople - dSLE_air  - LE_AC_Heat ...
                        + AC_onDehum*10^6*(qbin - q_RHspCooling); % (W/m2) building ground area
else
    YBuildInt(1)	=	TemperatureB(1) - 273.15;
    YBuildInt(2)	=	TemperatureB(2) - 273.15;
    YBuildInt(3)	=	TemperatureB(3) - 273.15;
    YBuildInt(4)	=	TemperatureB(4) - 273.15;
    YBuildInt(5)	=	TemperatureB(5) - 273.15;
    YBuildInt(6)    =   TemperatureB(6) - 273.15;
    YBuildInt(7)    =   TemperatureB(7) - 273.15;
    YBuildInt(8)    =   1000*TemperatureB(8) - 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Waste heat into the canyon, total anthropogenic heat input due to cooling
% and heating, and building energy demand
if AC_on==1
    % Positive if AC is on and added to the canyon as sensible and latent heat source
    % If the moisture removed from the air is fully condensed, we need to
    % account for tis additional removal of latent heat in the sensible
    % heat flux
    % The sensible and latent heat flux from the AC does not only account
    % for heat brought in through ventilation, but also heat that is
    % internally released by the walls
    WasteHeat.SensibleFromAC_Can    = (H_AC_Heat + LE_AC_Heat + (H_AC_Heat+LE_AC_Heat)/COP_AC)*Wroof/Wcan; % (W/m^2) canyon ground
    WasteHeat.LatentFromAC_Can      = LE_AC_Heat*Wroof/Wcan; %  W/m^2 canyon ground
    WasteHeat.WaterFromAC_Can       = LE_AC_Heat*Wroof/Wcan; %  W/m^2 canyon ground, water that is condensed and removed as runoff in sewer system 
    
    WasteHeat.SensibleFromHeat_Can  = 0; % W/m^2 canyon ground
    WasteHeat.LatentFromHeat_Can    = 0; % W/m^2 canyon ground
    
    WasteHeat.SensibleFromVent_Can  = -Hvent*Wroof/Wcan; % W/m^2 canyon ground: it is negative for the canyon as hot air is leaving for the cooler indoor air
    WasteHeat.LatentFromVent_Can    = -LEvent*Wroof/Wcan; % W/m^2 canyon ground

    % Anthropogenic heat added to the urban area accounts for
    % the internal sources in the building of sensible and latent heat due
    % to equipment and people and the additional heat added to the system
    % because of the coefficient of performance AC not being infinite (some additional energy is
    % needed to move the heat from one place to the other. 
    WasteHeat.TotAnthInput_URB      =   ((H_AC_Heat+LE_AC_Heat)/COP_AC + Hequip + Hpeople  + LEequip + LEpeople)*Wroof/(Wcan+Wroof); % W/m^2 urban area

    EnergyUse.EnergyForAC       = dth*(H_AC_Heat+LE_AC_Heat)*Wroof/COP_AC; % W*h
    EnergyUse.EnergyForAC_H     = dth*(H_AC_Heat)*Wroof/COP_AC; % W*h
    EnergyUse.EnergyForAC_LE    = dth*(LE_AC_Heat)*Wroof/COP_AC; % W*h
    EnergyUse.EnergyForHeating  = 0; % W*h

elseif  Heat_on==1
    WasteHeat.SensibleFromAC_Can    = 0; % W/m^2 canyon ground
    WasteHeat.LatentFromAC_Can      = 0; % W/m^2 canyon ground
    WasteHeat.WaterFromAC_Can       = 0; % W/m^2 canyon ground
    WasteHeat.SensibleFromHeat_Can  = 0; % W/m^2 canyon ground, all waste heat from heating is added to the canyon through ventilation and the conductive heat fluxes through the walls
    WasteHeat.LatentFromHeat_Can    = 0; % W/m^2 canyon ground

    WasteHeat.SensibleFromVent_Can  = -Hvent*Wroof/Wcan; % W/m^2 canyon ground, positive -> heat added to the canyon
    WasteHeat.LatentFromVent_Can    = -LEvent*Wroof/Wcan; % W/m^2 canyon ground, positive -> heat added to the canyon

    % Anthropogenic heat added to the urban area accounts for
    % the internal sources in the building of sensible and latent heat due
    % to equipment and people and the additional heat added to the system needed for heating
    % COP for eating is not considered in the waseheat (adding heat is the
    % intended purpose of the heating system and this is usually done in
    % the building)
    WasteHeat.TotAnthInput_URB      = (-H_AC_Heat-LE_AC_Heat + Hequip + Hpeople + LEequip + LEpeople)*Wroof/(Wcan+Wroof);

    EnergyUse.EnergyForAC       = 0; % W*h
    EnergyUse.EnergyForAC_H     = 0; % W*h
    EnergyUse.EnergyForAC_LE    = 0; % W*h
    EnergyUse.EnergyForHeating  = dth*(-H_AC_Heat-LE_AC_Heat)*Wroof/COP_Heat; % W*h
    
else
    WasteHeat.SensibleFromAC_Can    = 0; % W/m^2 canyon ground
    WasteHeat.LatentFromAC_Can      = 0; % W/m^2 canyon ground
    WasteHeat.WaterFromAC_Can       = 0; % W/m^2 canyon ground
    WasteHeat.SensibleFromHeat_Can  = 0; % W/m^2 canyon ground
    WasteHeat.LatentFromHeat_Can    = 0; % W/m^2 canyon ground

    WasteHeat.SensibleFromVent_Can  = -Hvent*Wroof/Wcan; % W/m^2 canyon ground, positive -> heat added to the canyon
    WasteHeat.LatentFromVent_Can    = -LEvent*Wroof/Wcan; % W/m^2 canyon ground, positive -> heat added to the canyon

    WasteHeat.TotAnthInput_URB      = (Hequip + Hpeople + LEequip + LEpeople)*Wroof/(Wcan+Wroof); % W/m^2 canyon ground

    EnergyUse.EnergyForAC       = 0; % W*h, COP for heating and cooling could be different
    EnergyUse.EnergyForAC_H     = 0; % W*h
    EnergyUse.EnergyForAC_LE    = 0; % W*h
    EnergyUse.EnergyForHeating  = 0; % W*h
end

% Apply fraction of Air conditioned rooms for calculation of energy use and
% for waste heat calculation
WasteHeat.SensibleFromAC_Can    = HVACSchedule.AirConRoomFraction.*WasteHeat.SensibleFromAC_Can;
WasteHeat.LatentFromAC_Can      = HVACSchedule.AirConRoomFraction.*WasteHeat.LatentFromAC_Can;
WasteHeat.WaterFromAC_Can       = HVACSchedule.AirConRoomFraction.*WasteHeat.WaterFromAC_Can;
WasteHeat.SensibleFromHeat_Can  = HVACSchedule.AirConRoomFraction.*WasteHeat.SensibleFromHeat_Can;
WasteHeat.LatentFromHeat_Can    = HVACSchedule.AirConRoomFraction.*WasteHeat.LatentFromHeat_Can;
WasteHeat.SensibleFromVent_Can  = HVACSchedule.AirConRoomFraction.*WasteHeat.SensibleFromVent_Can;
WasteHeat.LatentFromVent_Can    = HVACSchedule.AirConRoomFraction.*WasteHeat.LatentFromVent_Can;
WasteHeat.TotAnthInput_URB      = HVACSchedule.AirConRoomFraction.*WasteHeat.TotAnthInput_URB;
EnergyUse.EnergyForAC           = HVACSchedule.AirConRoomFraction.*EnergyUse.EnergyForAC;
EnergyUse.EnergyForAC_H         = HVACSchedule.AirConRoomFraction.*EnergyUse.EnergyForAC_H;
EnergyUse.EnergyForAC_LE        = HVACSchedule.AirConRoomFraction.*EnergyUse.EnergyForAC_LE;
EnergyUse.EnergyForHeating      = HVACSchedule.AirConRoomFraction.*EnergyUse.EnergyForHeating;



% ----------- In case of no BEM
if BEM_on ~=1
    WasteHeat.SensibleFromAC_Can    = 0; % W/m^2 canyon ground
    WasteHeat.LatentFromAC_Can      = 0; % W/m^2 canyon ground
    WasteHeat.WaterFromAC_Can       = 0; % W/m^2 canyon ground
    WasteHeat.SensibleFromHeat_Can  = 0; % W/m^2 canyon ground
    WasteHeat.LatentFromHeat_Can    = 0; % W/m^2 canyon ground
    WasteHeat.SensibleFromVent_Can  = 0; % W/m^2 canyon ground, positive -> heat added to the canyon
    WasteHeat.LatentFromVent_Can    = 0; % W/m^2 canyon ground, positive -> heat added to the canyon
    WasteHeat.TotAnthInput_URB      = 0; % W/m^2 canyon ground
    EnergyUse.EnergyForAC       = 0; % W*h, COP for heating and cooling could be different
    EnergyUse.EnergyForAC_H     = 0; % W*h
    EnergyUse.EnergyForAC_LE    = 0; % W*h
    EnergyUse.EnergyForHeating  = 0; % W*h  
end



