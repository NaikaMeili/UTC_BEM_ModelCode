function[EnergyFluxUrban,EnergyFluxCan,EnergyFluxRoof]=PlanAreaEnergyBalanceCalculation(ViewFactor,MeteoDataRaw,...
    SWRin,SWRout,SWRabs,LWRin,LWRout,LWRabs,LEflux,Hflux,Gflux,...
    geometry_Out,FractionsGround_Out,PropOpticalRoof_Out,Anthropo,Figure)


% Incoming radiation
%--------------------------------------------------------------------------
SWRin_atm			=	(MeteoDataRaw.SAB1_in + MeteoDataRaw.SAB2_in + MeteoDataRaw.SAD1_in + MeteoDataRaw.SAD2_in);
SWRinDir_atm		=	(MeteoDataRaw.SAB1_in + MeteoDataRaw.SAB2_in);
SWRinDiff_atm		=	(MeteoDataRaw.SAD1_in + MeteoDataRaw.SAD2_in);

LWRin_atm			=	MeteoDataRaw.LWR_in;

% View factors
%--------------------------------------------------------------------------
F_gs_T=ViewFactor.F_gs_T; F_gt_T=ViewFactor.F_gt_T; F_gw_T=ViewFactor.F_gw_T;
F_ww_T=ViewFactor.F_ww_T; F_wt_T=ViewFactor.F_wt_T; F_wg_T=ViewFactor.F_wg_T;
F_ws_T=ViewFactor.F_ws_T; F_sg_T=ViewFactor.F_sg_T; F_sw_T=ViewFactor.F_sw_T;
F_st_T=ViewFactor.F_st_T; F_tg_T=ViewFactor.F_tg_T; F_tw_T=ViewFactor.F_tw_T;
F_ts_T=ViewFactor.F_ts_T; F_tt_T=ViewFactor.F_tt_T;

% Normalize surface area
%--------------------------------------------------------------------------
A_s     =   geometry_Out.wcanyon; 
A_g     =   geometry_Out.wcanyon; 
A_w     =   geometry_Out.hcanyon;
A_t     =   2*2*pi*geometry_Out.radius_tree;	% There are 2 trees. Hence, the area of tree is twice a circle
fgveg   =   FractionsGround_Out.fveg; 
fgbare  =   FractionsGround_Out.fbare; 
fgimp   =   FractionsGround_Out.fimp;

% Rescaling tree absorbed radiation
SWRabs.SWRabsTree = SWRabs.SWRabsTree.*4.*geometry_Out.radius_tree./(4.*geometry_Out.radius_tree.*pi);
LWRabs.LWRabsTree = LWRabs.LWRabsTree.*4.*geometry_Out.radius_tree./(4.*geometry_Out.radius_tree.*pi);

% Shortwave canyon components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CanSWRin_SurfArea	=	SWRin.SWRinGroundVeg.*fgveg.*A_g./A_g + SWRin.SWRinGroundBare.*fgbare*A_g./A_g + SWRin.SWRinGroundImp.*fgimp.*A_g./A_g + ...
						SWRin.SWRinWallSun.*A_w./A_g + SWRin.SWRinWallShade.*A_w./A_g + SWRin.SWRinTree.*A_t./A_g;

CanSWRabs_SurfArea	=	SWRabs.SWRabsGroundVeg.*fgveg.*A_g./A_g + SWRabs.SWRabsGroundBare*fgbare.*A_g./A_g + SWRabs.SWRabsGroundImp*fgimp.*A_g./A_g + ...
						SWRabs.SWRabsWallSun.*A_w./A_g + SWRabs.SWRabsWallShade.*A_w./A_g + SWRabs.SWRabsTree.*A_t./A_g;
      
CanSWRout_SurfArea	=	SWRout.SWRoutGroundVeg.*fgveg.*A_g/A_s+SWRout.SWRoutGroundBare*fgbare.*A_g./A_s+SWRout.SWRoutGroundImp.*fgimp.*A_g./A_s+...
						SWRout.SWRoutWallSun.*A_w./A_s+SWRout.SWRoutWallShade.*A_w./A_s+SWRout.SWRoutTree.*A_t./A_s;
	  
CanSWRout_Ref_to_Atm=	SWRout.SWRoutGroundVeg.*F_sg_T.*fgveg + SWRout.SWRoutGroundBare.*F_sg_T.*fgbare + SWRout.SWRoutGroundImp.*F_sg_T.*fgimp + ...
						SWRout.SWRoutWallSun.*F_sw_T + SWRout.SWRoutWallShade.*F_sw_T + SWRout.SWRoutTree.*F_st_T;
					
CanSWR_EB_SurfArea			=	CanSWRin_SurfArea - CanSWRabs_SurfArea - CanSWRout_SurfArea;
CanSWR_EB_PlanAreaCanyon	=	SWRin_atm - CanSWRabs_SurfArea - CanSWRout_Ref_to_Atm;


% Shortwave urban components
%--------------------------------------------------------------------------
% Calculate total urban including roofs
UrbanSWRin_SurfArea     =   geometry_Out.wcanyon_norm.*CanSWRin_SurfArea + geometry_Out.wroof_norm*SWRin.SWRinTotalRoof;
UrbanSWRabs_SurfArea    =   geometry_Out.wcanyon_norm.*CanSWRabs_SurfArea + geometry_Out.wroof_norm*SWRabs.SWRabsTotalRoof;
UrbanSWRout_SurfArea    =   geometry_Out.wcanyon_norm.*CanSWRout_SurfArea + geometry_Out.wroof_norm*SWRout.SWRoutTotalRoof;
UrbanSWRout_Ref_to_Atm  =   geometry_Out.wcanyon_norm.*CanSWRout_Ref_to_Atm + geometry_Out.wroof_norm*SWRout.SWRoutTotalRoof;

UrbanSWR_EB_SurfArea			=	UrbanSWRin_SurfArea - UrbanSWRabs_SurfArea - UrbanSWRout_SurfArea;
UrbanSWR_EB_PlanAreaUrban	=	SWRin_atm - UrbanSWRabs_SurfArea - UrbanSWRout_Ref_to_Atm;

% Canyon albedo calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UrbanAlbedo	=	UrbanSWRout_Ref_to_Atm./SWRin_atm;
CanAlbedo	=	CanSWRout_Ref_to_Atm./SWRin_atm;
RoofAlbedo	=	PropOpticalRoof_Out.albedo;


% Longwave canyon components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CanLWRin_SurfArea	=	LWRin.LWRinGroundVeg.*fgveg.*A_g./A_g + LWRin.LWRinGroundBare.*fgbare*A_g./A_g + LWRin.LWRinGroundImp.*fgimp.*A_g./A_g + ...
						LWRin.LWRinWallSun.*A_w./A_g + LWRin.LWRinWallShade.*A_w./A_g + LWRin.LWRinTree.*A_t./A_g;

CanLWRabs_SurfArea	=	LWRabs.LWRabsGroundVeg.*fgveg.*A_g./A_g + LWRabs.LWRabsGroundBare*fgbare.*A_g./A_g + LWRabs.LWRabsGroundImp*fgimp.*A_g./A_g + ...
						LWRabs.LWRabsWallSun.*A_w./A_g + LWRabs.LWRabsWallShade.*A_w./A_g + LWRabs.LWRabsTree.*A_t./A_g;
      
CanLWRout_SurfArea	=	LWRout.LWRoutGroundVeg.*fgveg.*A_g/A_s+LWRout.LWRoutGroundBare*fgbare.*A_g./A_s+LWRout.LWRoutGroundImp.*fgimp.*A_g./A_s+...
						LWRout.LWRoutWallSun.*A_w./A_s+LWRout.LWRoutWallShade.*A_w./A_s+LWRout.LWRoutTree.*A_t./A_s;
	  
CanLWRout_Ref_to_Atm=	LWRout.LWRoutGroundVeg.*F_sg_T.*fgveg + LWRout.LWRoutGroundBare.*F_sg_T.*fgbare + LWRout.LWRoutGroundImp.*F_sg_T.*fgimp + ...
						LWRout.LWRoutWallSun.*F_sw_T + LWRout.LWRoutWallShade.*F_sw_T + LWRout.LWRoutTree.*F_st_T;
					
CanLWR_EB_SurfArea			=	CanLWRin_SurfArea - CanLWRabs_SurfArea - CanLWRout_SurfArea;
CanLWR_EB_PlanAreaCanyon	=	LWRin_atm - CanLWRabs_SurfArea - CanLWRout_Ref_to_Atm;


% Longwave urban components
%--------------------------------------------------------------------------
% Calculate total urban including roofs
UrbanLWRin_SurfArea     =   geometry_Out.wcanyon_norm.*CanLWRin_SurfArea + geometry_Out.wroof_norm*LWRin.LWRinTotalRoof;
UrbanLWRabs_SurfArea    =   geometry_Out.wcanyon_norm.*CanLWRabs_SurfArea + geometry_Out.wroof_norm*LWRabs.LWRabsTotalRoof;
UrbanLWRout_SurfArea    =   geometry_Out.wcanyon_norm.*CanLWRout_SurfArea + geometry_Out.wroof_norm*LWRout.LWRoutTotalRoof;
UrbanLWRout_Ref_to_Atm  =   geometry_Out.wcanyon_norm.*CanLWRout_Ref_to_Atm + geometry_Out.wroof_norm*LWRout.LWRoutTotalRoof;

UrbanLWR_EB_SurfArea			=	UrbanLWRin_SurfArea - UrbanLWRabs_SurfArea - UrbanLWRout_SurfArea;
UrbanLWR_EB_PlanAreaUrban	=	LWRin_atm - UrbanLWRabs_SurfArea - UrbanLWRout_Ref_to_Atm;


% Urban energy balance components
%--------------------------------------------------------------------------
Qanth		=	Anthropo.Qf_canyon.*geometry_Out.wcanyon_norm + Anthropo.Qf_roof.*geometry_Out.wcanyon_norm; % Anthropogenic heat flux (upward)


EBtot = SWRabs.SWRabsTotalUrban + LWRabs.LWRabsTotalUrban -...
    LEflux.LEfluxUrban - Hflux.HfluxUrban - Gflux.G1Urban + Qanth;

% Critical Energy Balance components out
%--------------------------------------------------------------------------
% Urban
EnergyFluxUrban.SWRin_PlanArea          = SWRin_atm;
EnergyFluxUrban.SWRin_SurfArea          = UrbanSWRin_SurfArea;
EnergyFluxUrban.SWRabs_SurfArea         = UrbanSWRabs_SurfArea;
EnergyFluxUrban.SWRout_SurfArea         = UrbanSWRout_SurfArea;
EnergyFluxUrban.SWRout_Ref_to_Atm       = UrbanSWRout_Ref_to_Atm;
EnergyFluxUrban.SWREB_SurfArea          = UrbanSWR_EB_SurfArea;
EnergyFluxUrban.SWREB_PlanAreaUrban     = UrbanSWR_EB_PlanAreaUrban;

EnergyFluxUrban.LWRin_PlanArea          = LWRin_atm;
EnergyFluxUrban.LWRin_SurfArea          = UrbanLWRin_SurfArea;
EnergyFluxUrban.LWRabs_SurfArea         = UrbanLWRabs_SurfArea;
EnergyFluxUrban.LWRout_SurfArea         = UrbanLWRout_SurfArea;
EnergyFluxUrban.LWRout_Ref_to_Atm       = UrbanLWRout_Ref_to_Atm;
EnergyFluxUrban.LWREB_SurfArea          = UrbanLWR_EB_SurfArea;
EnergyFluxUrban.LWREB_PlanAreaUrban     = UrbanLWR_EB_PlanAreaUrban;

EnergyFluxUrban.UrbanAlbedo  = UrbanAlbedo;

EnergyFluxUrban.SWRabs	= SWRabs.SWRabsTotalUrban;
EnergyFluxUrban.LWRabs	= LWRabs.LWRabsTotalUrban;
EnergyFluxUrban.LEflux	= LEflux.LEfluxUrban;
EnergyFluxUrban.Hflux	= Hflux.HfluxUrban;
EnergyFluxUrban.Gflux	= Gflux.G1Urban;
EnergyFluxUrban.Qanth   = Qanth;
EnergyFluxUrban.EB      = EBtot;

% Canyon
EnergyFluxCan.SWRin_SurfArea          = CanSWRin_SurfArea;
EnergyFluxCan.SWRabs_SurfArea         = CanSWRabs_SurfArea;
EnergyFluxCan.SWRout_SurfArea         = CanSWRout_SurfArea;
EnergyFluxCan.SWRout_Ref_to_Atm       = CanSWRout_Ref_to_Atm;
EnergyFluxCan.SWREB_SurfArea          = CanSWR_EB_SurfArea;
EnergyFluxCan.SWREB_PlanAreaCanyon	= CanSWR_EB_PlanAreaCanyon;

EnergyFluxCan.LWRin_SurfArea          = CanLWRin_SurfArea;
EnergyFluxCan.LWRabs_SurfArea         = CanLWRabs_SurfArea;
EnergyFluxCan.LWRout_SurfArea         = CanLWRout_SurfArea;
EnergyFluxCan.LWRout_Ref_to_Atm       = CanLWRout_Ref_to_Atm;
EnergyFluxCan.LWREB_SurfArea          = CanLWR_EB_SurfArea;
EnergyFluxCan.LWREB_PlanAreaCanyon  = CanLWR_EB_PlanAreaCanyon;

EnergyFluxCan.CanAlbedo  = CanAlbedo;

EnergyFluxCan.SWRabs	= SWRabs.SWRabsTotalCanyon;
EnergyFluxCan.LWRabs	= LWRabs.LWRabsTotalCanyon;
EnergyFluxCan.LEflux	= LEflux.LEfluxCanyon;
EnergyFluxCan.Hflux     = Hflux.HfluxCanyon;
EnergyFluxCan.Gflux     = Gflux.G1Canyon;
EnergyFluxCan.Qanth     = Anthropo.Qf_canyon;
EnergyFluxCan.EB        = EnergyFluxCan.SWRabs + EnergyFluxCan.LWRabs - EnergyFluxCan.LEflux - ...
    EnergyFluxCan.Hflux - EnergyFluxCan.Gflux + EnergyFluxCan.Qanth;

% Roof
EnergyFluxRoof.SWRin	= SWRin.SWRinTotalRoof;
EnergyFluxRoof.SWRout	= SWRout.SWRoutTotalRoof;
EnergyFluxRoof.LWRin	= LWRin.LWRinTotalRoof;
EnergyFluxRoof.LWRout	= LWRout.LWRoutTotalRoof;

EnergyFluxRoof.RoofAlbedo  = RoofAlbedo;

EnergyFluxRoof.SWRabs	= SWRabs.SWRabsTotalRoof;
EnergyFluxRoof.LWRabs	= LWRabs.LWRabsTotalRoof;
EnergyFluxRoof.LEflux	= LEflux.LEfluxRoof;
EnergyFluxRoof.Hflux    = Hflux.HfluxRoof;
EnergyFluxRoof.Gflux    = Gflux.G1Roof;
EnergyFluxRoof.Qanth    = Anthropo.Qf_roof;
EnergyFluxRoof.EB       = EnergyFluxRoof.SWRabs + EnergyFluxRoof.LWRabs - EnergyFluxRoof.LEflux - ...
    EnergyFluxRoof.Hflux - EnergyFluxRoof.Gflux + EnergyFluxRoof.Qanth;

EnergyFluxRoof.SWREB	= SWRin.SWRinTotalRoof - SWRabs.SWRabsTotalRoof - SWRout.SWRoutTotalRoof;
EnergyFluxRoof.LWREB    = LWRin.LWRinTotalRoof - LWRabs.LWRabsTotalRoof - LWRout.LWRoutTotalRoof;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot graphs
if Figure==1
    
TTUrban = table(hour(MeteoDataRaw.Date),month(MeteoDataRaw.Date),...
            EnergyFluxUrban.SWRin_PlanArea,EnergyFluxUrban.SWRabs_SurfArea,EnergyFluxUrban.SWRout_Ref_to_Atm,...
            EnergyFluxUrban.LWRin_PlanArea,EnergyFluxUrban.LWRabs_SurfArea,EnergyFluxUrban.LWRout_Ref_to_Atm,...
            EnergyFluxUrban.LEflux,EnergyFluxUrban.Hflux,EnergyFluxUrban.Gflux,EnergyFluxUrban.Qanth,...
            EnergyFluxUrban.EB,EnergyFluxUrban.UrbanAlbedo,...
            EnergyFluxUrban.Hflux./EnergyFluxUrban.LEflux);

TTUrban.Properties.VariableNames = {'Hour','Month','SWRin','SWRabs','SWRout','LWRin','LWRabs','LWRout',...
    'LE','H','G','Qanth','EB','Albedo','BowenRatio'};

TTUrbanDiurnal = varfun(@nanmean,TTUrban,'GroupingVariables','Hour');
TTUrbanSeasonal = varfun(@nanmean,TTUrban,'GroupingVariables','Month');

TTUrbanDiurnalMedian = varfun(@median,TTUrban,'GroupingVariables','Hour');
TTUrbanSeasonalMedian = varfun(@median,TTUrban,'GroupingVariables','Month');


% Energy budget
f1 = figure;
set(f1, 'Units','centimeters','Position', [1 1 19 7])

t = tiledlayout(1,3);
t.Padding = 'compact'; %t.TileSpacing = 'compact';

nexttile
plot(MeteoDataRaw.Date,EnergyFluxUrban.EB,'k','DisplayName','EB')
xlabel('time'); ylabel('EB W/m^{2}'); title('Time series');
%legend

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_EB,'k','LineWidth',1.5,'DisplayName','EB')
xlim([0 23]); xlabel('hour'); ylabel('EB W/m^{2}'); title('Diurnal');
%legend

nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_EB,'k','LineWidth',1.5,'DisplayName','EB')
xlim([1 12]); xlabel('Month'); ylabel('EB W/m^{2}'); title('Seasonal');
%legend
sgtitle('Energy budget closure')

% Albedo und Bowen ratio
f1 = figure;
set(f1, 'Units','centimeters','Position', [1 1 14 10])

t = tiledlayout(2,2);
t.Padding = 'compact'; %t.TileSpacing = 'compact';

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_Albedo,'k','LineWidth',1.5,'DisplayName','albedo')
xlim([0 23]); xlabel('hour'); ylabel('Albedo (-)'); title('Albedo = SWR_{out}/SWR_{in}');
subtitle('Diurnal');
%legend

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnalMedian.median_BowenRatio,'k','LineWidth',1.5,'DisplayName','Bowen ratio')
xlim([0 23]); xlabel('hour'); ylabel('Bowen Ratio (-)'); title('Bowen ratio = H/LE, median');
subtitle('Diurnal');
%legend

nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_Albedo,'k','LineWidth',1.5,'DisplayName','albedo')
xlim([1 12]); xlabel('month'); ylabel('Albedo (-)'); subtitle('Seasonal');
%legend

nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonalMedian.median_BowenRatio,'k','LineWidth',1.5,'DisplayName','Bowen ratio')
xlim([1 12]); xlabel('month'); ylabel('Bowen Ratio (-)'); subtitle('Seasonal');
%legend

% Energy Fluxes
f1 = figure;
set(f1, 'Units','centimeters','Position', [1 1 19 14])

t = tiledlayout(2,3);
t.Padding = 'compact'; %t.TileSpacing = 'compact';

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_SWRin,'k','LineWidth',1.5,'DisplayName','SWR_{in}')
hold on
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_SWRabs,'r','LineWidth',1.5,'DisplayName','SWR_{abs}')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_SWRout,'b','LineWidth',1.5,'DisplayName','SWR_{out}')
xlim([0 23]); xlabel('hour'); ylabel('SWR W/m^{2}'); title('SWR'); subtitle('Diurnal');
legend('Location','southoutside','NumColumns',2)

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_LWRin,'k','LineWidth',1.5,'DisplayName','LWR_{in}')
hold on
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_LWRabs,'r','LineWidth',1.5,'DisplayName','LWR_{abs}')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_LWRout,'b','LineWidth',1.5,'DisplayName','LWR_{out}')
xlim([0 23]); xlabel('hour'); ylabel('LWR W/m^{2}'); title('LWR'); subtitle('Diurnal');
legend('Location','southoutside','NumColumns',2)

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_SWRabs,'k','LineWidth',1.5,'DisplayName','SWR_{abs}')
hold on
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_H,'r','LineWidth',1.5,'DisplayName','H')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_LE,'g','LineWidth',1.5,'DisplayName','LE')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_G,'b','LineWidth',1.5,'DisplayName','G')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_LWRabs,'m','LineWidth',1.5,'DisplayName','LWR_{abs}')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_Qanth,'k:','LineWidth',1.5,'DisplayName','Q_{f}')
xlim([0 23]); xlabel('hour'); ylabel('EB W/m^{2}'); title('Energy fluxes'); subtitle('Diurnal');
legend('Location','southoutside','NumColumns',2)

nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_SWRin,'k','LineWidth',1.5,'DisplayName','SWR_{in}')
hold on
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_SWRabs,'r','LineWidth',1.5,'DisplayName','SWR_{abs}')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_SWRout,'b','LineWidth',1.5,'DisplayName','SWR_{out}')
xlim([1 12]); xlabel('month'); ylabel('SWR W/m^{2}'); subtitle('Seasonal'); %title('SWR'); 
%legend

nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_LWRin,'k','LineWidth',1.5,'DisplayName','LWR_{in}')
hold on
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_LWRabs,'r','LineWidth',1.5,'DisplayName','LWR_{abs}')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_LWRout,'b','LineWidth',1.5,'DisplayName','LWR_{out}')
xlim([1 12]); xlabel('month'); ylabel('LWR W/m^{2}');  subtitle('Seasonal'); %title('LWR');
%legend

nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_SWRabs,'k','LineWidth',1.5,'DisplayName','SWR_{abs}')
hold on
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_H,'r','LineWidth',1.5,'DisplayName','H')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_LE,'g','LineWidth',1.5,'DisplayName','LE')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_G,'b','LineWidth',1.5,'DisplayName','G')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_LWRabs,'m','LineWidth',1.5,'DisplayName','LWR_{abs}')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_Qanth,'k:','LineWidth',1.5,'DisplayName','Q_{f}')
xlim([1 12]); xlabel('month'); ylabel('EB W/m^{2}'); subtitle('Seasonal'); %title('EB');
%legend
sgtitle('Energy fluxes')

end
