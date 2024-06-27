function[WaterFluxRoof,WaterFluxCan,WaterFluxBuild,WaterFluxUrban]=WaterBalanceComponents(MeteoDataRaw,...
    Runon,Leakage,LEflux,dVwater_dt,OwaterInitial,Owater,dInt_dt,Int,Anthropo,...
    LEbuildInt,BEMWasteHeat,ParHVAC,...
    ParSoilRoof_Out,ParSoilGround_Out,ParCalculation_Out,FractionsRoof_Out,...
    FractionsGround_Out,geometry_Out,Gemeotry_m_Out,Figure,ittmTot)


for ittm=1:ittmTot

% Calculation parameters
L_heat			=	1000.*(2501.3 - 2.361.*(MeteoDataRaw.T_atm-273.15));		% Latent heat vaporization/condensation [J/kg]
dts				=	ParCalculation_Out.dts;		% time step of calculation [s]
dth				=	ParCalculation_Out.dth;

% Water balance components
%--------------------------------------------------------------------------
% Precipitation
%--------------------------------------------------------------------------
RainRoof	=   MeteoDataRaw.rain; %[mm/time step]
RainCan     =   MeteoDataRaw.rain; %[mm/time step]
RainUrb     =   MeteoDataRaw.rain; %[mm/time step]

% Irrigation at the surface
%--------------------------------------------------------------------------
IrrSurfRoof	=   FractionsRoof_Out(ittm).fveg.*Anthropo.Qf_roof(:,1,ittm); %[mm/time step]
IrrSurfCan	=   FractionsGround_Out(ittm).fveg.*Anthropo.Waterf_canyonVeg(:,1,ittm) + ...
                FractionsGround_Out(ittm).fbare.*Anthropo.Waterf_canyonBare(:,1,ittm); %[mm/time step]
IrrSurfUrb	=   geometry_Out(ittm).wcanyon_norm.*IrrSurfCan + geometry_Out(ittm).wroof_norm.*IrrSurfRoof; %[mm/time step]

% Runoff leaving the system
%--------------------------------------------------------------------------
RunoffRoof	=	Runon.RunoffRoofTot(:,1,ittm);    %[mm/time step] 
RunoffCan	=	Runon.RunoffGroundTot(:,1,ittm);  %[mm/time step] 
RunoffUrb	=	Runon.RunoffUrban(:,1,ittm);      %[mm/time step]

% Leakage at bottom of soil column, leaving the system
%--------------------------------------------------------------------------
LeakageRoof	=	dth.*Leakage.LkRoof(:,1,ittm);	% [mm/time step] Subsurface runoff (out of gridcell)
LeakageCan	=	dth.*Leakage.LkGround(:,1,ittm);	% [mm/time step] Subsurface runoff (out of gridcell)
LeakageUrb	=	dth.*Leakage.LkUrban(:,1,ittm);	% [mm/time step] Subsurface runoff (out of gridcell)

% Evapotranspiration, latent heat
%--------------------------------------------------------------------------
ETRoof      =   LEflux.LEfluxRoof(:,1,ittm)./L_heat.*dts;% Total evapotranspiration (upward)
ETCan       =   LEflux.LEfluxCanyon(:,1,ittm)./L_heat.*dts;% Total evapotranspiration (upward)
ETUrb       =   LEflux.LEfluxUrban(:,1,ittm)./L_heat.*dts;% Total evapotranspiration (upward)

% Change in latent heat stored in the air
%--------------------------------------------------------------------------
dS_ET_dtBuild   =   LEbuildInt.dSLE_air(:,1,ittm)./L_heat.*dts;
dS_ET_dtCan     =   LEflux.dS_LE_air(:,1,ittm)./L_heat.*dts;
dS_ET_dtUrb     =   geometry_Out(ittm).wroof_norm.*dS_ET_dtBuild + geometry_Out(ittm).wcanyon_norm.*dS_ET_dtCan;

% Urban evapotranspiration fluxes according to their source: 1) Interception
% and ponding, 2) evaporation from soil, 3) transpiration, 4) Latent heat
% sources in building interior
%--------------------------------------------------------------------------
ETEvapoIntUrb   =   (geometry_Out(ittm).wroof_norm.*(FractionsRoof_Out(ittm).fimp.*LEflux.LEfluxRoofImp(:,1,ittm) + ...
                    FractionsRoof_Out(ittm).fveg.*(LEflux.LEfluxRoofVegInt(:,1,ittm) + LEflux.LEfluxRoofVegPond(:,1,ittm))) + ...
                    geometry_Out(ittm).wcanyon_norm.*(FractionsGround_Out(ittm).fimp.*LEflux.LEfluxGroundImp(:,1,ittm) + ...
                    FractionsGround_Out(ittm).fbare.*LEflux.LEfluxGroundBarePond(:,1,ittm) + ...
                    FractionsGround_Out(ittm).fveg.*(LEflux.LEfluxGroundVegInt(:,1,ittm) + LEflux.LEfluxGroundVegPond(:,1,ittm)) + ...
                     4.*geometry_Out(ittm).radius_tree.*LEflux.LEfluxTreeInt(:,1,ittm)))./L_heat.*dts;

ETEvapoSoilUrb  =   (geometry_Out(ittm).wroof_norm.*FractionsRoof_Out(ittm).fveg.*LEflux.LEfluxRoofVegSoil(:,1,ittm) + ...
                    geometry_Out(ittm).wcanyon_norm.*(FractionsGround_Out(ittm).fbare.* LEflux.LEfluxGroundBareSoil(:,1,ittm) + ...
                    FractionsGround_Out(ittm).fveg.*LEflux.LEfluxGroundVegSoil(:,1,ittm)))./L_heat.*dts;
                
ETTranspUrb     =   (geometry_Out(ittm).wroof_norm.*FractionsRoof_Out(ittm).fveg.*LEflux.LTEfluxRoofVeg(:,1,ittm) + ...
                    geometry_Out(ittm).wcanyon_norm.*(FractionsGround_Out(ittm).fveg.*LEflux.LTEfluxGroundVeg(:,1,ittm) + ...
                    4.*geometry_Out(ittm).radius_tree.*LEflux.LTEfluxTree(:,1,ittm)))./L_heat.*dts;

ET_buildAnth    =   geometry_Out(ittm).wroof_norm.*(LEbuildInt.LEpeople(:,1,ittm) + LEbuildInt.LEequip(:,1,ittm))./L_heat.*dts;

                                
% Change in soil moisture
%--------------------------------------------------------------------------
dVdtRoof    =   FractionsRoof_Out(ittm).fveg.*dVwater_dt.dVRoofSoilVeg_dt(:,1,ittm);
dVdtCan     =   dVwater_dt.dVGroundSoilTot_dt(:,1,ittm);
dVdtUrb     =	geometry_Out(ittm).wcanyon_norm.*dVdtCan + geometry_Out(ittm).wroof_norm.*dVdtRoof;

[dVdtRoofCalc,dVdtCanCalc,dVdtUrbCalc]=soil_functions.PostCalculateSoilMoistureChange(...
    OwaterInitial,Owater,ParSoilRoof_Out,ParSoilGround_Out,FractionsRoof_Out,FractionsGround_Out,geometry_Out,ittm);

% Irrigation within soil (due to fixed soil moisture)
%--------------------------------------------------------------------------
IrrSoilRoof	=   dVdtRoofCalc - dVdtRoof;
IrrSoilCan  =   dVdtCanCalc - dVdtCan;
IrrSoilUrb  =   dVdtUrbCalc - dVdtUrb;

% Change in intercepted water
%--------------------------------------------------------------------------
% On plant canopy
dIdtPlantRoof   =   FractionsRoof_Out(ittm).fveg.*dInt_dt.dInt_dtRoofVegPlant(:,1,ittm);
dIdtPlantCan    =   FractionsGround_Out(ittm).fveg.*dInt_dt.dInt_dtGroundVegPlant(:,1,ittm)+...
                    4.*geometry_Out(ittm).radius_tree.*dInt_dt.dInt_dtTree(:,1,ittm);
dIdtPlantUrb    =   geometry_Out(ittm).wcanyon_norm.*dIdtPlantCan + geometry_Out(ittm).wroof_norm.*dIdtPlantRoof;
% On ground/surface
dIdtGroundRoof  =   FractionsRoof_Out(ittm).fveg.*dInt_dt.dInt_dtRoofVegGround(:,1,ittm) +... 
                    FractionsRoof_Out(ittm).fimp.*dInt_dt.dInt_dtRoofImp(:,1,ittm);
dIdtGroundCan   =   FractionsGround_Out(ittm).fveg.*dInt_dt.dInt_dtGroundVegGround(:,1,ittm) +...
                    FractionsGround_Out(ittm).fbare.*dInt_dt.dInt_dtGroundBare(:,1,ittm) +...
                    FractionsGround_Out(ittm).fimp.*dInt_dt.dInt_dtGroundImp(:,1,ittm);
dIdtGroundUrb	=   geometry_Out(ittm).wcanyon_norm.*dIdtGroundCan + geometry_Out(ittm).wroof_norm.*dIdtGroundRoof;
% Due to runon
dRun_dtRoof		=	Runon.RunonRoofTot(:,1,ittm) - [0; Runon.RunonRoofTot(1:end-1,1,ittm)];
dRun_dtCan		=	Runon.RunonGroundTot(:,1,ittm) - [0; Runon.RunonGroundTot(1:end-1,1,ittm)];
dRun_dtUrb		=	Runon.RunonUrban(:,1,ittm) - [0; Runon.RunonUrban(1:end-1,1,ittm)];
% Total
dIdtRoof	=   dIdtPlantRoof + dIdtGroundRoof + dRun_dtRoof;
dIdtCan     =   dIdtPlantCan + dIdtGroundCan + dRun_dtCan;
dIdtUrb     =   dIdtPlantUrb + dIdtGroundUrb + dRun_dtUrb;

% Surface water storage (SurfStor)
%--------------------------------------------------------------------------
IntRoof     =   Int.IntRooftot(:,1,ittm) + Runon.RunonRoofTot(:,1,ittm);
IntCan      =   FractionsGround_Out(ittm).fimp.*Int.IntGroundImp(:,1,ittm) + ...
                FractionsGround_Out(ittm).fbare.*Int.IntGroundBare(:,1,ittm) + ...
                FractionsGround_Out(ittm).fveg.*(Int.IntGroundVegPlant(:,1,ittm) + Int.IntGroundVegGround(:,1,ittm)) +...
                4.*geometry_Out(ittm).radius_tree.*Int.IntTree(:,1,ittm) + Runon.RunonGroundTot(:,1,ittm);
IntUrb      =   geometry_Out(ittm).wcanyon_norm.*IntCan + geometry_Out(ittm).wroof_norm.*IntRoof;

IntUrb_tm1  =   [0; IntUrb(1:end-1)];
dIdtUrbCalc =   IntUrb - IntUrb_tm1;


% Anthropogenic latent heat fluxes due to sources in building and HVAC
%--------------------------------------------------------------------------
% Anthropogenic sources in building interior: people & equipment
ET_buildAnthRoof    =   0;
ET_buildAnthCan     =   0;
ET_buildAnthBuild   =  (LEbuildInt.LEpeople(:,1,ittm) + LEbuildInt.LEequip(:,1,ittm))./L_heat.*dts;
ET_buildAnthUrb     =   geometry_Out(ittm).wroof_norm.*(LEbuildInt.LEpeople(:,1,ittm) + LEbuildInt.LEequip(:,1,ittm))./L_heat.*dts;

% Latent heat removed from air due to ventilation
ET_VentRoof         =   0;
ET_VentCan          =   BEMWasteHeat.LatentFromVent_Can(:,1,ittm)./L_heat.*dts;
ET_VentBuild        =   LEbuildInt.LEvent(:,1,ittm)./L_heat.*dts;
ET_VentUrb          =   0;

% ET exchange due to AC
ET_ACRoof   =   0;
ET_ACCan    =   BEMWasteHeat.LatentFromAC_Can(:,1,ittm)./L_heat.*dts;
ET_ACBuild  =   BEMWasteHeat.LatentFromAC_Can(:,1,ittm).*(Gemeotry_m_Out(ittm).Width_canyon./Gemeotry_m_Out(ittm).Width_roof)./L_heat.*dts;
ET_ACUrb    =   0;

% ET exchange due to HVAC between indoor and outdoor air
ET_HVACexchRoof   =   0;
ET_HVACexchCan    =   ET_VentCan + ET_ACCan;
ET_HVACexchBuild  =   ET_VentBuild - ET_ACBuild;
ET_HVACexchUrb    =   0;

% Water removed from building interior due to condensation during air
% conditioning
ET_WasteWaterACBuild    =   BEMWasteHeat.WaterFromAC_Can(:,1,ittm).*(Gemeotry_m_Out(ittm).Width_canyon./Gemeotry_m_Out(ittm).Width_roof)./L_heat.*dts;

ET_WasteWaterACUrb      =   BEMWasteHeat.WaterFromAC_Can(:,1,ittm).*Gemeotry_m_Out(ittm).Width_canyon./(Gemeotry_m_Out(ittm).Width_canyon+Gemeotry_m_Out(ittm).Width_roof)./L_heat.*dts;



% Water balance
%--------------------------------------------------------------------------
WBRoof  =   RainRoof + IrrSurfRoof + IrrSoilRoof - RunoffRoof - LeakageRoof - ETRoof - dVdtRoofCalc - dIdtRoof;
WBCan   =   RainCan + IrrSurfCan + IrrSoilCan - RunoffCan - LeakageCan - ETCan - dVdtCanCalc - dIdtCan - dS_ET_dtCan; % ET_HVACexchCan is already includd in ETCan;
WBBuild =   ET_buildAnthBuild + ET_HVACexchBuild - ET_WasteWaterACBuild - dS_ET_dtBuild;
WBUrb   =   RainUrb + IrrSurfUrb + IrrSoilUrb + ET_buildAnthUrb - RunoffUrb - LeakageUrb - ETUrb - dVdtUrbCalc - dIdtUrb - dS_ET_dtUrb; % Moisture removed from the air is already accounted for in ETUrb


% figure
% tiledlayout(2,2)
% nexttile; plot(WBRoof); title('Roof');
% nexttile; plot(WBCan); title('Canyon');
% nexttile; plot(WBUrb); title('Urban');
% nexttile; plot(WBBuild); title('Building');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Incoming rainfall
WaterFluxRoof.Rain(:,1,ittm)	=   RainRoof; %[mm/time step]
WaterFluxCan.Rain(:,1,ittm)	=   RainCan; %[mm/time step]
WaterFluxUrban.Rain(:,1,ittm)	=   RainUrb; %[mm/time step]

% Runoff leaving the system
WaterFluxRoof.Runoff(:,1,ittm)	=   RunoffRoof; %[mm/time step]
WaterFluxCan.Runoff(:,1,ittm)     =   RunoffCan; %[mm/time step]
WaterFluxUrban.Runoff(:,1,ittm)	=   RunoffUrb; %[mm/time step]

% Leakage at the bottom of the soil column -> water leaving the system
WaterFluxRoof.Leakage(:,1,ittm)	=   LeakageRoof; %[mm/time step]
WaterFluxCan.Leakage(:,1,ittm)	=   LeakageCan; %[mm/time step]
WaterFluxUrban.Leakage(:,1,ittm)	=   LeakageUrb; %[mm/time step]

% Evapotranspiration flux
WaterFluxRoof.ET(:,1,ittm)	=   ETRoof; %[mm/time step]
WaterFluxCan.ET(:,1,ittm)     =   ETCan; %[mm/time step]
WaterFluxUrban.ET(:,1,ittm)	=   ETUrb; %[mm/time step]
% Evapotranspiration based on different sources of the moisture flux
WaterFluxUrban.ETEvaporationFromSurface(:,1,ittm)	=   ETEvapoIntUrb; %[mm/time step]
WaterFluxUrban.ETEvaporationFromSoil(:,1,ittm)	=   ETEvapoSoilUrb; %[mm/time step]
WaterFluxUrban.ETTranspiration(:,1,ittm)          =   ETTranspUrb; %[mm/time step]
WaterFluxUrban.ETAnthSourceBuildInt(:,1,ittm)     =   ET_buildAnth; %[mm/time step]

% Change in moisture storage in the air
WaterFluxBuild.dS_ET_dt(:,1,ittm) =   dS_ET_dtBuild; %[mm/time step]
WaterFluxCan.dS_ET_dt(:,1,ittm)	=   dS_ET_dtCan; %[mm/time step]
WaterFluxUrban.dS_ET_dt(:,1,ittm)	=   dS_ET_dtUrb; %[mm/time step]

% Exchange of water vapor between indoor and outdoor air
WaterFluxRoof.ET_HVACexch(:,1,ittm)   =   ET_HVACexchRoof; %[mm/time step]
WaterFluxCan.ET_HVACexch(:,1,ittm)    =   ET_HVACexchCan; %[mm/time step]
WaterFluxBuild.ET_HVACexch(:,1,ittm)  =   ET_HVACexchBuild; %[mm/time step]
WaterFluxUrban.ET_HVACexch(:,1,ittm)  =   ET_HVACexchUrb; %[mm/time step]

% Anthropogenic moisture flux due to people/equipment in buildings 
WaterFluxRoof.AnthBuildInt(:,1,ittm)  =   0; %[mm/time step]
WaterFluxCan.AnthBuildInt(:,1,ittm)   =   0; %[mm/time step]
WaterFluxBuild.AnthBuildInt(:,1,ittm) =   ET_buildAnthBuild; %[mm/time step]
WaterFluxUrban.AnthBuildInt(:,1,ittm)	=   ET_buildAnthUrb; %[mm/time step]

% Anthropogenic removal of moisture from air due to condensation during
% air-conditioning process
WaterFluxRoof.WasteWaterAC(:,1,ittm)  =   0; %[mm/time step]
WaterFluxCan.WasteWaterAC(:,1,ittm)   =   0; %[mm/time step]
WaterFluxBuild.WasteWaterAC(:,1,ittm) =   ET_WasteWaterACBuild; %[mm/time step]
WaterFluxUrban.WasteWaterAC(:,1,ittm)	=   ET_WasteWaterACUrb; %[mm/time step]

% Change in water storage in the soil due to soil moisture change
WaterFluxRoof.dVdt(:,1,ittm)	=   dVdtRoofCalc; %[mm/time step]
WaterFluxCan.dVdt(:,1,ittm)	=   dVdtCanCalc; %[mm/time step]
WaterFluxUrban.dVdt(:,1,ittm)	=   dVdtUrbCalc; %[mm/time step]

% Chang in water storage due to interception change
WaterFluxRoof.dIdt(:,1,ittm)	=   dIdtRoof; %[mm/time step]
WaterFluxCan.dIdt(:,1,ittm)	=   dIdtCan; %[mm/time step]
WaterFluxUrban.dIdt(:,1,ittm)	=   dIdtUrb; %[mm/time step]

% Irrigation flux applied at the surface of the soil
WaterFluxRoof.IrrSurf(:,1,ittm)	=   IrrSurfRoof; %[mm/time step]
WaterFluxCan.IrrSurf(:,1,ittm)	=   IrrSurfCan; %[mm/time step]
WaterFluxUrban.IrrSurf(:,1,ittm)	=   IrrSurfUrb; %[mm/time step]

% Irrigation occuring due to fixed soil moisture in some soil layers
WaterFluxRoof.IrrSoil(:,1,ittm)	=   IrrSoilRoof; %[mm/time step]
WaterFluxCan.IrrSoil(:,1,ittm)	=   IrrSoilCan; %[mm/time step]
WaterFluxUrban.IrrSoil(:,1,ittm)	=   IrrSoilUrb; %[mm/time step]

% Total irritation flux (includes surface irrigation and fixed soil
% moisture in certain soil layers)
WaterFluxRoof.IrrTot(:,1,ittm)	=   IrrSoilRoof + IrrSurfRoof; %[mm/time step]
WaterFluxCan.IrrTot(:,1,ittm)     =   IrrSoilCan + IrrSurfCan; %[mm/time step]
WaterFluxUrban.IrrTot(:,1,ittm)	=   IrrSoilUrb + IrrSurfUrb; %[mm/time step]

% Intercepted or ponding water
WaterFluxRoof.Int(:,1,ittm)	=   IntRoof; %[mm]
WaterFluxCan.Int(:,1,ittm)	=   IntCan; %[mm]
WaterFluxUrban.Int(:,1,ittm)	=   IntUrb; %[mm]

% Water balance
WaterFluxRoof.WB(:,1,ittm)	=   WBRoof; %[mm/time step]
WaterFluxBuild.WB(:,1,ittm)   =   WBBuild; %[mm/time step]
WaterFluxCan.WB(:,1,ittm)     =   WBCan; %[mm/time step]
WaterFluxUrban.WB(:,1,ittm)	=   WBUrb; %[mm/time step]

% Water Balance second test
%--------------------------------------------------------------------------
WBRoof2     =   WaterFluxRoof.Rain(:,1,ittm) + WaterFluxRoof.IrrTot(:,1,ittm) - WaterFluxRoof.Runoff(:,1,ittm) - WaterFluxRoof.Leakage(:,1,ittm)...
                - WaterFluxRoof.ET(:,1,ittm) - WaterFluxRoof.dIdt(:,1,ittm) - WaterFluxRoof.dVdt(:,1,ittm);

WBBuild2    =   WaterFluxBuild.AnthBuildInt(:,1,ittm) + WaterFluxBuild.ET_HVACexch(:,1,ittm) - WaterFluxBuild.WasteWaterAC(:,1,ittm) - WaterFluxBuild.dS_ET_dt(:,1,ittm);

WBCan2      =   WaterFluxCan.Rain(:,1,ittm) + WaterFluxCan.IrrTot(:,1,ittm) - WaterFluxCan.Runoff(:,1,ittm) - WaterFluxCan.Leakage(:,1,ittm)...
                -  WaterFluxCan.ET(:,1,ittm) - WaterFluxCan.dIdt(:,1,ittm) - WaterFluxCan.dVdt(:,1,ittm) - WaterFluxCan.dS_ET_dt(:,1,ittm); % ET_HVACexchCan is already includd in ETCan;

WBUrb2      =   WaterFluxUrban.Rain(:,1,ittm) + WaterFluxUrban.IrrTot(:,1,ittm) - WaterFluxUrban.Runoff(:,1,ittm) - WaterFluxUrban.Leakage(:,1,ittm)...
                -WaterFluxUrban.ET(:,1,ittm) - WaterFluxUrban.dIdt(:,1,ittm) - WaterFluxUrban.dVdt(:,1,ittm) - WaterFluxUrban.dS_ET_dt(:,1,ittm)...
                + WaterFluxUrban.AnthBuildInt(:,1,ittm); % Moisture removed from the air due AC condensation is already accounted for in ETUrb

% figure
% tiledlayout(2,2)
% nexttile; plot(WBRoof2(2:end)); title('Roof');
% nexttile; plot(WBCan2(2:end)); title('Canyon');
% nexttile; plot(WBUrb2(2:end)); title('Urban');
% nexttile; plot(WBBuild2(2:end)); title('Building');
% sgtitle(['Water balance, ittm = ' num2str(ittm)])

% ET Test
% ETbalance2 = WaterFluxUrban.ET(:,1,ittm) - (WaterFluxUrban.ETEvaporationFromSurface(:,1,ittm) + ...
%             WaterFluxUrban.ETEvaporationFromSoil(:,1,ittm) + WaterFluxUrban.ETTranspiration(:,1,ittm)...
%             + WaterFluxUrban.AnthBuildInt(:,1,ittm) - WaterFluxUrban.dS_ET_dt(:,1,ittm));
% 
% figure
% plot(ETbalance2); title('ET balance');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Figure==1

for ittm=1:ittmTot    
% Plot graphs
TTUrban = table(hour(MeteoDataRaw.Date),month(MeteoDataRaw.Date),...
            WaterFluxUrban.Rain(:,1,ittm),WaterFluxUrban.Runoff(:,1,ittm),WaterFluxUrban.Leakage(:,1,ittm),...
            WaterFluxUrban.ET(:,1,ittm),WaterFluxUrban.ETEvaporationFromSurface(:,1,ittm),...
            WaterFluxUrban.ETEvaporationFromSoil(:,1,ittm),WaterFluxUrban.ETTranspiration(:,1,ittm),...
            WaterFluxUrban.dVdt(:,1,ittm),WaterFluxUrban.dIdt(:,1,ittm),WaterFluxUrban.IrrSurf(:,1,ittm),WaterFluxUrban.IrrSoil(:,1,ittm),...
            WaterFluxUrban.IrrTot(:,1,ittm),WaterFluxUrban.Int(:,1,ittm),WaterFluxUrban.WB(:,1,ittm),...
            WaterFluxUrban.dS_ET_dt(:,1,ittm),WaterFluxUrban.AnthBuildInt(:,1,ittm),WaterFluxUrban.WasteWaterAC(:,1,ittm));
        
TTUrban.Properties.VariableNames = {'Hour','Month','Rain','Runoff','Leakage',...
    'ET','ETEvaporationFromSurface','ETEvaporationFromSoil','ETTranspiration',...
    'dVdt','dIdt','IrrSurf','IrrSoil','IrrTot','Int','WB',...
    'dSETdt','QanthBuild','WasteWaterAC'};

TTUrbanDiurnal = varfun(@nanmean,TTUrban,'GroupingVariables','Hour');
TTUrbanSeasonal = varfun(@nanmean,TTUrban,'GroupingVariables','Month');


% Water budget WB
f1 = figure;
set(f1, 'Units','centimeters','Position', [1 1 19 7])

t = tiledlayout(1,3);
t.Padding = 'compact'; %t.TileSpacing = 'compact';

nexttile
plot(MeteoDataRaw.Date,WaterFluxUrban.WB(:,1,ittm),'k','DisplayName','WB')
xlabel('time'); ylabel('WB mm/time step'); title('Time series');
%legend

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_WB,'k','LineWidth',1.5,'DisplayName','WB')
xlim([0 23]); xlabel('hour'); ylabel('WB mm/time step'); title('Diurnal');
%legend

nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_WB,'k','LineWidth',1.5,'DisplayName','WB')
xlim([1 12]); xlabel('Month'); ylabel('WB mm/time step'); title('Seasonal');
%legend
sgtitle(['Water budget closure, ittm = ' num2str(ittm)])

% Evapotranspiration and interception
f1 = figure;
set(f1, 'Units','centimeters','Position', [1 1 14 7])

t = tiledlayout(1,2);
t.Padding = 'compact'; %t.TileSpacing = 'compact';

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_ET,'k','LineWidth',1.5,'DisplayName','ET_{tot}')
hold on
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_ETEvaporationFromSurface,'r','LineWidth',1.5,'DisplayName','E_{surface}')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_ETEvaporationFromSoil,'b','LineWidth',1.5,'DisplayName','E_{soil}')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_ETTranspiration,'g','LineWidth',1.5,'DisplayName','T_{vegetation}')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_QanthBuild,'y','LineWidth',1.5,'DisplayName','Q_{LE,anth,building}')
plot(TTUrbanDiurnal.Hour,-TTUrbanDiurnal.nanmean_WasteWaterAC,'m','LineWidth',1.5,'DisplayName','Q_{AC,waste water}')
xlim([0 23]); xlabel('hour'); ylabel('ET mm/time step'); subtitle('Diurnal');

nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_ET,'k','LineWidth',1.5,'DisplayName','ET_{tot}')
hold on
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_ETEvaporationFromSurface,'r','LineWidth',1.5,'DisplayName','E_{surface}')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_ETEvaporationFromSoil,'b','LineWidth',1.5,'DisplayName','E_{soil}')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_ETTranspiration,'g','LineWidth',1.5,'DisplayName','T_{vegetation}')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_QanthBuild,'y','LineWidth',1.5,'DisplayName','Q_{LE,anth,building}')
plot(TTUrbanSeasonal.Month,-TTUrbanSeasonal.nanmean_WasteWaterAC,'m','LineWidth',1.5,'DisplayName','Q_{AC,waste water}')
xlim([1 12]); xlabel('Month'); ylabel('ET mm/time step'); subtitle('Seasonal');
legend('Location','NorthEastOutside')
sgtitle(['Evapotranspiration flux partitioning, ittm = ' num2str(ittm)])


% Water fluxes
f1 = figure;
set(f1, 'Units','centimeters','Position', [1 1 19 8])

t = tiledlayout(1,2);
t.Padding = 'compact'; %t.TileSpacing = 'compact';

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_Rain,'k','LineWidth',1.5,'DisplayName','Rain')
hold on
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_IrrTot,'k--','LineWidth',1.5,'DisplayName','Irrigation')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_Runoff,'b','LineWidth',1.5,'DisplayName','Runoff')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_ET,'g','LineWidth',1.5,'DisplayName','ET')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_dVdt,'r','LineWidth',1.5,'DisplayName','dVdt')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_dIdt,'m','LineWidth',1.5,'DisplayName','dIdt')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_Leakage,'c','LineWidth',1.5,'DisplayName','Leakage')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_QanthBuild,'y','LineWidth',1.5,'DisplayName','Q_{LE,anth,build}')
xlim([0 23]); xlabel('hour'); ylabel('Water fluxes mm/time step'); subtitle('Diurnal');

nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_Rain,'k','LineWidth',1.5,'DisplayName','Rain')
hold on
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_IrrTot,'k--','LineWidth',1.5,'DisplayName','Irrigation')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_Runoff,'b','LineWidth',1.5,'DisplayName','Runoff')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_ET,'g','LineWidth',1.5,'DisplayName','ET')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_dVdt,'r','LineWidth',1.5,'DisplayName','dVdt')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_dIdt,'m','LineWidth',1.5,'DisplayName','dIdt')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_Leakage,'c','LineWidth',1.5,'DisplayName','Leakage')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_QanthBuild,'y','LineWidth',1.5,'DisplayName','Q_{LE,anth,build}')
xlim([1 12]); xlabel('Month'); ylabel('Water fluxes mm/time step'); subtitle('Seasonal');
legend('Location','NorthEastOutside')
sgtitle(['Water fluxes, ittm = ' num2str(ittm)])

end
end


