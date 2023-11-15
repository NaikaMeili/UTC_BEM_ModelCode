function[HbinWallSun,HbinWallshd,HBinRoof,HBinGround,HbinIntMass,HbinWindow]=...
    SensibleHeatFluxBuildingInterior(Tbin,Tinwallsun,Tinwallshd,Tceiling,Tinground,Tintmass,Twindow)


% Turbulent sensible heat flux from walls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The convective heat transfer coefficient according to Bueno et al. 2012,
% Oleson et a. 2019
hvert   =   3.076; % (W/m^2K) vertical (walls) 
hhorRed =   0.948; % (W/m^2K) for a horizontal surface with reduced convection (floor surface with Tsi <Tin and ceiling surface with Tsi > Tin)
hhorEnh =   4.040; % (W/m^2K) for a horizontal surface with enhanced convection (floor surface with Tsi > Tin and ceiling surface with Tsi < Tin)

% Sensible heat flux from wall and roof into building air
HbinWallSun     =   hvert.*(Tinwallsun-Tbin); % (W/m^2), sensible heat from sunlit wall to building interior 
HbinWallshd     =   hvert.*(Tinwallshd-Tbin); % (W/m^2), sensible heat from sunlit wall to building interior
HbinWindow      =   hvert.*(Twindow-Tbin); % (W/m^2), sensible heat from window to building interior
HbinIntMass     =   hvert.*(Tintmass-Tbin); % (W/m^2), sensible heat from internal mass (e.g. internal walls) to building interior
HBinRoof        =   ((Tceiling>Tbin).*hhorRed + (Tceiling<Tbin).*hhorEnh).*(Tceiling-Tbin); % (W/m^2), sensible heat from impervious roof to building interior
HBinGround      =   ((Tinground>Tbin).*hhorEnh + (Tinground<Tbin).*hhorRed).*(Tinground-Tbin); % (W/m^2), sensible heat from ground to building interior (assume to be zero)





