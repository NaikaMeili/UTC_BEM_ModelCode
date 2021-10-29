function[WaterFluxRoof,WaterFluxCan,WaterFluxUrban]=WaterBalanceComponents(MeteoDataRaw,...
    Runon,Leakage,LEflux,dVwater_dt,OwaterInitial,Owater,dInt_dt,Int,Anthropo,...
    ParSoil,ParCalculation_Out,FractionsRoof_Out,FractionsGround_Out,geometry_Out,Figure)

% Calculation parameters
L_heat			=	1000.*(2501.3 - 2.361.*(MeteoDataRaw.T_atm-273.15));		% Latent heat vaporization/condensation [J/kg]
dts				=	ParCalculation_Out.dts;		% time step of calculation [s]
dth				=	ParCalculation_Out.dth;

% Water balance components
%--------------------------------------------------------------------------
% Precipitation
RainRoof	=   MeteoDataRaw.rain; %[mm/time step]
RainCan     =   MeteoDataRaw.rain; %[mm/time step]
RainUrb     =   MeteoDataRaw.rain; %[mm/time step]

% Irrigation at the surface
IrrSurfRoof	=   FractionsRoof_Out.fveg.*Anthropo.Qf_roof; %[mm/time step]
IrrSurfCan	=   FractionsGround_Out.fveg.*Anthropo.Waterf_canyonVeg + ...
                FractionsGround_Out.fbare.*Anthropo.Waterf_canyonBare; %[mm/time step]
IrrSurfUrb	=   geometry_Out.wcanyon_norm.*IrrSurfCan + geometry_Out.wroof_norm.*IrrSurfRoof; %[mm/time step]

% Runoff leaving the system
RunoffRoof	=	Runon.RunoffRoofTot;    %[mm/time step] 
RunoffCan	=	Runon.RunoffGroundTot;  %[mm/time step] 
RunoffUrb	=	Runon.RunoffUrban;      %[mm/time step]

% Leakage at bottom of soil column, leaving the system
LeakageRoof	=	dth.*Leakage.LkRoof;	% [mm/time step] Subsurface runoff (out of gridcell)
LeakageCan	=	dth.*Leakage.LkGround;	% [mm/time step] Subsurface runoff (out of gridcell)
LeakageUrb	=	dth.*Leakage.LkUrban;	% [mm/time step] Subsurface runoff (out of gridcell)

% Evapotranspiration, latent heat
ETRoof      =   LEflux.LEfluxRoof./L_heat.*dts;% Total evapotranspiration (upward)
ETCan       =   LEflux.LEfluxCanyon./L_heat.*dts;% Total evapotranspiration (upward)
ETUrb       =   LEflux.LEfluxUrban./L_heat.*dts;% Total evapotranspiration (upward)

% Further ET fluxes, according to source
ETEvapoIntUrb   =   (geometry_Out.wroof_norm.*(FractionsRoof_Out.fimp.*LEflux.LEfluxRoofImp + ...
                    FractionsRoof_Out.fveg.*(LEflux.LEfluxRoofVegInt + LEflux.LEfluxRoofVegPond)) + ...
                    geometry_Out.wcanyon_norm.*(FractionsGround_Out.fimp.*LEflux.LEfluxGroundImp + ...
                    FractionsGround_Out.fbare.*LEflux.LEfluxGroundBarePond + ...
                    FractionsGround_Out.fveg.*(LEflux.LEfluxGroundVegInt + LEflux.LEfluxGroundVegPond) + ...
                     4.*geometry_Out.radius_tree.*LEflux.LEfluxTreeInt))./L_heat.*dts;

ETEvapoSoilUrb  =   (geometry_Out.wroof_norm.*FractionsRoof_Out.fveg.*LEflux.LEfluxRoofVegSoil + ...
                    geometry_Out.wcanyon_norm.*(FractionsGround_Out.fbare.* LEflux.LEfluxGroundBareSoil + ...
                    FractionsGround_Out.fveg.*LEflux.LEfluxGroundVegSoil))./L_heat.*dts;
                
ETTranspUrb     =   (geometry_Out.wroof_norm.*FractionsRoof_Out.fveg.*LEflux.LTEfluxRoofVeg + ...
                    geometry_Out.wcanyon_norm.*(FractionsGround_Out.fveg.*LEflux.LTEfluxGroundVeg + ...
                    4.*geometry_Out.radius_tree.*LEflux.LTEfluxTree))./L_heat.*dts;

ETTest = ETUrb - (ETEvapoIntUrb + ETEvapoSoilUrb + ETTranspUrb);
                
                
% Change in soil moisture
dVdtRoof    =   FractionsRoof_Out.fveg.*dVwater_dt.dVRoofSoilVeg_dt;
dVdtCan     =   dVwater_dt.dVGroundSoilTot_dt;
dVdtUrb     =	geometry_Out.wcanyon_norm.*dVdtCan + geometry_Out.wroof_norm.*dVdtRoof;

[dVdtRoofCalc,dVdtCanCalc,dVdtUrbCalc]=soil_functions.PostCalculateSoilMoistureChange(...
    OwaterInitial,Owater,ParSoil,FractionsRoof_Out,FractionsGround_Out,geometry_Out);

% Irrigation within soil (due to fixed soil moisture)
IrrSoilRoof	=   dVdtRoofCalc - dVdtRoof;
IrrSoilCan  =   dVdtCanCalc - dVdtCan;
IrrSoilUrb  =   dVdtUrbCalc - dVdtUrb;

% Change in intercepted water
% On plant canopy
dIdtPlantRoof   =   FractionsRoof_Out.fveg.*dInt_dt.dInt_dtRoofVegPlant;
dIdtPlantCan    =   FractionsGround_Out.fveg.*dInt_dt.dInt_dtGroundVegPlant+...
                    4.*geometry_Out.radius_tree.*dInt_dt.dInt_dtTree;
dIdtPlantUrb    =   geometry_Out.wcanyon_norm.*dIdtPlantCan + geometry_Out.wroof_norm.*dIdtPlantRoof;
% On ground/surface
dIdtGroundRoof  =   FractionsRoof_Out.fveg.*dInt_dt.dInt_dtRoofVegGround +... 
                    FractionsRoof_Out.fimp.*dInt_dt.dInt_dtRoofImp;
dIdtGroundCan   =   FractionsGround_Out.fveg.*dInt_dt.dInt_dtGroundVegGround +...
                    FractionsGround_Out.fbare.*dInt_dt.dInt_dtGroundBare +...
                    FractionsGround_Out.fimp.*dInt_dt.dInt_dtGroundImp;
dIdtGroundUrb	=   geometry_Out.wcanyon_norm.*dIdtGroundCan + geometry_Out.wroof_norm.*dIdtGroundRoof;
% Due to runon
dRun_dtRoof		=	Runon.RunonRoofTot - [0; Runon.RunonRoofTot(1:end-1)];
dRun_dtCan		=	Runon.RunonGroundTot - [0; Runon.RunonGroundTot(1:end-1)];
dRun_dtUrb		=	Runon.RunonUrban - [0; Runon.RunonUrban(1:end-1)];
% Total
dIdtRoof	=   dIdtPlantRoof + dIdtGroundRoof + dRun_dtRoof;
dIdtCan     =   dIdtPlantCan + dIdtGroundCan + dRun_dtCan;
dIdtUrb     =   dIdtPlantUrb + dIdtGroundUrb + dRun_dtUrb;

% surface water storage (SurfStor)
IntRoof     =   Int.IntRooftot + Runon.RunonRoofTot;
IntCan      =   FractionsGround_Out.fimp.*Int.IntGroundImp + ...
                FractionsGround_Out.fbare.*Int.IntGroundBare + ...
                FractionsGround_Out.fveg.*(Int.IntGroundVegPlant + Int.IntGroundVegGround) +...
                4.*geometry_Out.radius_tree.*Int.IntTree + Runon.RunonGroundTot;
IntUrb      =   geometry_Out.wcanyon_norm.*IntCan + geometry_Out.wroof_norm.*IntRoof;

IntUrb_tm1  =   [0; IntUrb(1:end-1)];
dIdtUrbCalc =   IntUrb - IntUrb_tm1;

% Water balance
WBRoof  =   RainRoof + IrrSurfRoof + IrrSoilRoof - RunoffRoof - LeakageRoof - ETRoof - dVdtRoofCalc - dIdtRoof;
WBCan   =   RainCan + IrrSurfCan + IrrSoilCan - RunoffCan - LeakageCan - ETCan - dVdtCanCalc - dIdtCan;
WBUrb   =   RainUrb + IrrSurfUrb + IrrSoilUrb - RunoffUrb - LeakageUrb - ETUrb - dVdtUrbCalc - dIdtUrb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WaterFluxRoof.Rain	=   RainRoof; %[mm/time step]
WaterFluxCan.Rain	=   RainCan; %[mm/time step]
WaterFluxUrban.Rain	=   RainUrb; %[mm/time step]

WaterFluxRoof.Runoff	=   RunoffRoof; %[mm/time step]
WaterFluxCan.Runoff     =   RunoffCan; %[mm/time step]
WaterFluxUrban.Runoff	=   RunoffUrb; %[mm/time step]

WaterFluxRoof.Leakage	=   LeakageRoof; %[mm/time step]
WaterFluxCan.Leakage	=   LeakageCan; %[mm/time step]
WaterFluxUrban.Leakage	=   LeakageUrb; %[mm/time step]

WaterFluxRoof.ET	=   ETRoof; %[mm/time step]
WaterFluxCan.ET     =   ETCan; %[mm/time step]
WaterFluxUrban.ET	=   ETUrb; %[mm/time step]
WaterFluxUrban.ETEvaporationFromSurface	=   ETEvapoIntUrb; %[mm/time step]
WaterFluxUrban.ETEvaporationFromSoil	=   ETEvapoSoilUrb; %[mm/time step]
WaterFluxUrban.ETTranspiration          =   ETTranspUrb; %[mm/time step]

WaterFluxRoof.dVdt	=   dVdtRoofCalc; %[mm/time step]
WaterFluxCan.dVdt	=   dVdtCanCalc; %[mm/time step]
WaterFluxUrban.dVdt	=   dVdtUrbCalc; %[mm/time step]

WaterFluxRoof.dIdt	=   dIdtRoof; %[mm/time step]
WaterFluxCan.dIdt	=   dIdtCan; %[mm/time step]
WaterFluxUrban.dIdt	=   dIdtUrb; %[mm/time step]

WaterFluxRoof.IrrSurf	=   IrrSurfRoof; %[mm/time step]
WaterFluxCan.IrrSurf	=   IrrSurfCan; %[mm/time step]
WaterFluxUrban.IrrSurf	=   IrrSurfUrb; %[mm/time step]

WaterFluxRoof.IrrSoil	=   IrrSoilRoof; %[mm/time step]
WaterFluxCan.IrrSoil	=   IrrSoilCan; %[mm/time step]
WaterFluxUrban.IrrSoil	=   IrrSoilUrb; %[mm/time step]

WaterFluxRoof.IrrTot	=   IrrSoilRoof + IrrSurfRoof; %[mm/time step]
WaterFluxCan.IrrTot     =   IrrSoilCan + IrrSurfCan; %[mm/time step]
WaterFluxUrban.IrrTot	=   IrrSoilUrb + IrrSurfUrb; %[mm/time step]

WaterFluxRoof.Int	=   IntRoof; %[mm]
WaterFluxCan.Int	=   IntCan; %[mm]
WaterFluxUrban.Int	=   IntUrb; %[mm]

WaterFluxRoof.WB	=   WBRoof; %[mm]
WaterFluxCan.WB     =   WBCan; %[mm]
WaterFluxUrban.WB	=   WBUrb; %[mm]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% subplot(1,3,1)
% plot(WBRoof,'DisplayName','WB_{roof}')
% subplot(1,3,2)
% plot(WBCan,'DisplayName','WB_{can}')
% subplot(1,3,3)
% plot(WBUrb,'DisplayName','WB_{urb}')

if Figure==1
    
% Plot graphs
TTUrban = table(hour(MeteoDataRaw.Date),month(MeteoDataRaw.Date),...
            WaterFluxUrban.Rain,WaterFluxUrban.Runoff,WaterFluxUrban.Leakage,...
            WaterFluxUrban.ET,WaterFluxUrban.ETEvaporationFromSurface,...
            WaterFluxUrban.ETEvaporationFromSoil,WaterFluxUrban.ETTranspiration,...
            WaterFluxUrban.dVdt,WaterFluxUrban.dIdt,WaterFluxUrban.IrrSurf,WaterFluxUrban.IrrSoil,...
            WaterFluxUrban.IrrTot,WaterFluxUrban.Int,WaterFluxUrban.WB);
        



TTUrban.Properties.VariableNames = {'Hour','Month','Rain','Runoff','Leakage',...
    'ET','ETEvaporationFromSurface','ETEvaporationFromSoil','ETTranspiration',...
    'dVdt','dIdt','IrrSurf','IrrSoil','IrrTot','Int','WB'};

TTUrbanDiurnal = varfun(@nanmean,TTUrban,'GroupingVariables','Hour');
TTUrbanSeasonal = varfun(@nanmean,TTUrban,'GroupingVariables','Month');


% Water budget WB
f1 = figure;
set(f1, 'Units','centimeters','Position', [1 1 19 7])

t = tiledlayout(1,3);
t.Padding = 'compact'; %t.TileSpacing = 'compact';

nexttile
plot(MeteoDataRaw.Date,WaterFluxUrban.WB,'k','DisplayName','WB')
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
sgtitle('Water budget closure')

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
xlim([0 23]); xlabel('hour'); ylabel('ET mm/time step'); subtitle('Diurnal');

nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_ET,'k','LineWidth',1.5,'DisplayName','ET_{tot}')
hold on
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_ETEvaporationFromSurface,'r','LineWidth',1.5,'DisplayName','E_{surface}')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_ETEvaporationFromSoil,'b','LineWidth',1.5,'DisplayName','E_{soil}')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_ETTranspiration,'g','LineWidth',1.5,'DisplayName','T_{vegetation}')
xlim([1 12]); xlabel('Month'); ylabel('ET mm/time step'); subtitle('Seasonal');
legend('Location','NorthEastOutside')
sgtitle('Evapotranspiration flux partitioning')

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
xlim([1 12]); xlabel('Month'); ylabel('Water fluxes mm/time step'); subtitle('Seasonal');
legend('Location','NorthEastOutside')
sgtitle('Water fluxes')


end



