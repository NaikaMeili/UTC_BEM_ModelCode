%%%%%%%%%% RUN TIME SERIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile('+data_functions', 'ForcingData_ZH2010.mat'));

% Decide if a varying LAI timeseries is provided [1] or not [0]
[LAI_TimeSeries]=data_functions.VaryingLAIInput(0,'LAI_Zurich_Area'); 


n			=	size(MeteoDataZH_h,1);% Calculation length, there is no need to change this
m			=	1;					% Length for sensitivity analysis
Name_Site	=	'ZH';	% Name for Data_UEHM_site
Name_SiteFD	=	'ZH';		% Name for UEHMForcingData
OPTION_RAY	=	1; % Load precalculated view factors [1], Recalculate view factors [0]

NameOutput	=	'ZH_2010';


%% Meteo data
% LWR_in [W/m2], SAB1_in [W/m2], SAB2_in [W/m2], SAD1_in [W/m2], SAD2_in [W/m2]
% T_atm	[K], windspeed_u[m/s, pressure_atm [Pa], rain [mm/h], rel_humidity [-]
MeteoDataZH_h.Windspeed(MeteoDataZH_h.Windspeed(1:n,:)==0) = 0.01;	% Wind speed cannot be 0 otherwise the resistance function fails


MeteoDataRaw	=	struct('LWR_in',MeteoDataZH_h.LWRin(1:n,:),'SAB1_in',MeteoDataZH_h.SAB1(1:n,:),...
					'SAB2_in',MeteoDataZH_h.SAB2(1:n,:),'SAD1_in',MeteoDataZH_h.SAD1(1:n,:),...
					'SAD2_in',MeteoDataZH_h.SAD2(1:n,:),'T_atm',MeteoDataZH_h.Tatm(1:n,:),...
					'windspeed_u',MeteoDataZH_h.Windspeed(1:n,:),'pressure_atm',MeteoDataZH_h.Pressure_Pa(1:n,:),...
					'rain',MeteoDataZH_h.Precipitation(1:n,:),'rel_humidity',MeteoDataZH_h.RelativeHumidity(1:n,:)./100,...
					'Date',MeteoDataZH_h.Time(1:n,:));

%% Calculation starts here. No need to change anything after this point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
clearvars -except MeteoDataRaw n m Name_Site OPTION_RAY NameOutput Name_SiteFD LAI_TimeSeries

SoilWPot.SoilPotWGroundTot_H =0;
SoilWPot.SoilPotWGroundVeg_L=0;

% Meteo data at time step 1 for initialization				
[~,MeteoData,HumidityAtm,~,~,~]=feval(strcat('data_functions.UEHMForcingData_',Name_SiteFD),MeteoDataRaw,1,SoilWPot);

% Soil parameters
[~,~,~,~,~,~,ParSoilRoof,ParSoilGround,~,~,~,~,~,~,~,~,~,ParVegRoof,ParVegGround,ParVegTree,~]=feval(strcat('data_functions.Data_UEHM_site_',Name_Site),MeteoData,1,LAI_TimeSeries);

ParSoil		=	struct('Roof',ParSoilRoof,'Ground',ParSoilGround);

[~,~,~,ParSoil.Roof.Osat,ParSoil.Roof.Ohy,~,~,~,~,~,ParSoil.Roof.O33,...
    ~,~,~,~,~,~,~,~,~,~,~,~,~,~,~]...
	=soil_functions.SoilParametersTotal(ParSoilRoof.Pcla,ParSoilRoof.Psan,ParSoilRoof.Porg,...
    ParSoilRoof.Kfc,ParSoilRoof.Phy,ParSoilRoof.SPAR,ParSoilRoof.Kbot,...
	ParVegRoof.CASE_ROOT,ParVegRoof.CASE_ROOT,ParVegRoof.ZR95,ParVegRoof.ZR95,...
    ParVegRoof.ZR50,ParVegRoof.ZR50,ParVegRoof.ZRmax,ParVegRoof.ZRmax,ParSoilRoof.Zs);

[~,~,~,ParSoil.Ground.Osat,ParSoil.Ground.Ohy,~,~,~,~,~,ParSoil.Ground.O33,...
    ~,~,~,~,~,~,~,~,~,~,~,~,~,~,~]...
	=soil_functions.SoilParametersTotal(ParSoilGround.Pcla,ParSoilGround.Psan,ParSoilGround.Porg,...
    ParSoilGround.Kfc,ParSoilGround.Phy,ParSoilGround.SPAR,ParSoilGround.Kbot,...
	ParVegTree.CASE_ROOT,ParVegGround.CASE_ROOT,ParVegTree.ZR95,ParVegGround.ZR95,...
    ParVegTree.ZR50,ParVegGround.ZR50,ParVegTree.ZRmax,ParVegGround.ZRmax,ParSoilGround.Zs);

ParSoil.Roof.dz		=	diff(ParSoil.Roof.Zs);	% [mm]  Thickness of the soil layers
ParSoil.Ground.dz	=	diff(ParSoil.Ground.Zs);% [mm]  Thickness of the soil layers

ParSoil.Roof.Osat   = unique(ParSoil.Roof.Osat);
ParSoil.Roof.Ohy    = unique(ParSoil.Roof.Ohy);
ParSoil.Roof.O33    = unique(ParSoil.Roof.O33);

ParSoil.Ground.Osat = unique(ParSoil.Ground.Osat);
ParSoil.Ground.Ohy  = unique(ParSoil.Ground.Ohy);
ParSoil.Ground.O33  = unique(ParSoil.Ground.O33);

%% Initializing vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% Temperature: TempVec
%--------------------------------------------------------------------------Psi_soil_gveg
% TRoofImp		=	Temperature roof impervious area [K]
% TRoofVeg		=	Temperature roof vegetated area [K]
% TRoofIntImp	=	Interior temperature roof impervious area [K]
% TRoofIntVeg	=	Interior temperature roof vegetated area [K]
% TGroundImp	=	Temperature ground impervious area [K]
% TGroundBare	=	Temperature ground bare area [K]
% TGroundVeg	=	Temperature ground vegetated area [K]
% TTree			=	Temperature tree canopy [K]
% TWallSun		=	Temperature sunlit area [K]
% TWallShade	=	Temperature shaded area [K]
% TWallIntSun	=	Interior temperature sunlit wall [K]
% TWallIntShade	=	Interior temperature shaded wall [K]
% TCanyon		=	Temperature canyon [K]
% Tatm			=	Temperature atmosphere(measured) [K]

TempVecNames	=	{'TRoofImp';'TRoofVeg';'TRoofIntImp';'TRoofIntVeg';...
					'TGroundImp';'TGroundBare';'TGroundVeg';'TTree';'TWallSun';...
					'TWallShade';'TWallIntSun';'TWallIntShade';'TCanyon';'Tatm'};

for i=1:size(TempVecNames,1)
	TempVec.(cell2mat(TempVecNames(i)))			=	zeros(n,1,m);
	TempVec.(cell2mat(TempVecNames(i)))(1,:,:)	=	MeteoDataRaw.T_atm(i);
end

TempVec.Tatm	=	repmat(MeteoDataRaw.T_atm(:,1),1,1,m); % Temperature atmosphere(measured)

% Dampening temperature: TempDamp
%--------------------------------------------------------------------------
% TGroundImp	=	Dampening temperature ground impervious area [K]
% TGroundBare	=	Dampening temperature ground bare area [K]
% TGroundVeg	=	Dampening temperature ground vegetated area [K]
% TTree			=	Dampening temperature tree canopy [K]

TempDampNames	=	{'TDampGroundImp';'TDampGroundBare';'TDampGroundVeg';'TDampTree'};

for i=1:size(TempDampNames,1)
	TempDamp.(cell2mat(TempDampNames(i)))		=	zeros(n,1,m);
	TempDamp.(cell2mat(TempDampNames(i)))(1,:,:)=	nanmean(MeteoDataRaw.T_atm);
end

%% Humidity: Humidity
%--------------------------------------------------------------------------
% CanyonRelative: Relative humidity at canyon calculation height (-)
% CanyonSpecific: Specific humidity at canyon calculation height (kg/kg)
% CanyonVapourPre: Vapour pressure at canyon calculation height (Pa)
% CanyonRelativeSat: Saturation relative humidity at canyon calculation height (-), is always 1 
% CanyonSpecificSat: Specific humidity at saturation at canyon calculation height (kg/kg)
% CanyonVapourPreSat: Saturation vapor pressure at canyon calculation height (Pa)
% AtmRelative: Relative humidity at atmospheric forcing height (-)
% AtmSpecific: Specific humidity at atmospheric forcing height (kg/kg)
% AtmVapourPre: Vapor pressure at atmospheric forcing height (Pa)
% AtmRelativeSat: Saturation relative humidity at atmospheric forcing height (-), is always 1 
% AtmSpecificSat: Specific humidity at saturation at atmospheric forcing height (kg/kg)
% AtmVapourPreSat: Saturation vapour pressure at atmospheric forcing height (Pa)

HumidityNames	=	{'CanyonRelative';'CanyonSpecific';'CanyonVapourPre';'CanyonRelativeSat';...
					'CanyonSpecificSat';'CanyonVapourPreSat';'AtmRelative';'AtmSpecific';'AtmVapourPre';...
					'AtmRelativeSat';'AtmSpecificSat';'AtmVapourPreSat'};

for i=1:size(HumidityNames,1)
	Humidity.(cell2mat(HumidityNames(i)))	=	zeros(n,1,m);
end
Humidity.CanyonSpecific(1,:,:)				=	HumidityAtm.AtmSpecific;

%% Energy fluxes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Shortwave radiation
%--------------------------------------------------------------------------
% Shortwave radiation absorbed: SWRabs
%--------------------------------------------------------------------------
% SWRabsRoofImp		=	Absorbed shortwave radiation roof impervious area [W/m2 horizontal roof area]
% SWRabsRoofVeg		=	Absorbed shortwave radiation roof vegetated area [W/m2 horizontal roof area]
% SWRabsGroundImp	=	Absorbed shortwave radiation ground impervious area [W/m2 horizontal ground area]
% SWRabsGroundBare	=	Absorbed shortwave radiation ground bare area [W/m2 horizontal ground area]
% SWRabsGroundVeg	=	Absorbed shortwave radiation ground vegetated area [W/m2 horizontal ground area]
% SWRabsTree		=	Absorbed shortwave radiation tree canopy [W/m2 horizontally projected tree area: 4*radius]
% SWRabsWallSun		=	Absorbed shortwave radiation sunlit area [W/m2 vertical wall area]
% SWRabsWallShade	=	Absorbed shortwave radiation shaded area [W/m2 vertical wall area]
% SWRabsTotalRoof	=	Total absorbed shortwave radiation by the roof area [W/m2 horizontal roof area]
% SWRabsTotalGround	=	Total absorbed shortwave radiation by the canyon ground area [W/m2]
% SWRabsTotalCanyon	=	Total absorbed shortwave radiation by all the canyon facets [W/m2]
% SWRabsTotalUrban	=	Total absorbed shortwave radiation by all the urban elements (roof plus canyon) [W/m2]
%--------------------------------------------------------------------------
% Incoming shortwave radiation: SWRin
%--------------------------------------------------------------------------
% SWRinRoofImp		=	Incoming shortwave radiation roof impervious area [W/m2 horizontal roof area]
% SWRinRoofVeg		=	Incoming shortwave radiation roof vegetated area [W/m2 horizontal roof area]
% SWRinGroundImp	=	Incoming shortwave radiation ground impervious area [W/m2 horizontal ground area]
% SWRinGroundBare	=	Incoming shortwave radiation ground bare area [W/m2 horizontal ground area]
% SWRinGroundVeg	=	Incoming shortwave radiation ground vegetated area [W/m2 horizontal ground area]
% SWRinTree			=	Incoming shortwave radiation tree canopy [W/m2 horizontally projected tree area: 4*radius]
% SWRinWallSun		=	Incoming shortwave radiation sunlit area [W/m2 vertical wall area]
% SWRinWallShade	=	Incoming shortwave radiation shaded area [W/m2 vertical wall area]
% SWRinTotalRoof	=	Total incoming shortwave radiation by the roof area [W/m2 horizontal roof area]
% SWRinTotalGround	=	Total incoming shortwave radiation by the canyon ground area [W/m2]
% SWRinTotalCanyon	=	Total incoming shortwave radiation by all the canyon facets [W/m2]
% SWRinTotalUrban	=	Total incoming shortwave radiation by all the urban elements (roof plus canyon) [W/m2]
%--------------------------------------------------------------------------
% Outgoing shortwave radiation: SWRout
%--------------------------------------------------------------------------
% SWRoutRoofImp		=	Outgoing shortwave radiation roof impervious area [W/m2 horizontal roof area]
% SWRoutRoofVeg		=	Outgoing shortwave radiation roof vegetated area [W/m2 horizontal roof area]
% SWRoutGroundImp	=	Outgoing shortwave radiation ground impervious area [W/m2 horizontal ground area]
% SWRoutGroundBare	=	Outgoing shortwave radiation ground bare area [W/m2 horizontal ground area]
% SWRoutGroundVeg	=	Outgoing shortwave radiation ground vegetated area [W/m2 horizontal ground area]
% SWRoutTree		=	Outgoing shortwave radiation tree canopy [W/m2 horizontally projected tree area: 4*radius]
% SWRoutWallSun		=	Outgoing shortwave radiation sunlit area [W/m2 vertical wall area]
% SWRoutWallShade	=	Outgoing shortwave radiation shaded area [W/m2 vertical wall area]
% SWRoutTotalRoof	=	Total outgoing shortwave radiation by the roof area [W/m2 horizontal roof area]
% SWRoutTotalGround	=	Total outgoing shortwave radiation by the canyon ground area [W/m2]
% SWRoutTotalCanyon	=	Total outgoing shortwave radiation by all the canyon facets [W/m2]
% SWRoutTotalUrban	=	Total outgoing shortwave radiation by all the urban elements (roof plus canyon) [W/m2]
%--------------------------------------------------------------------------
% Shortwave radiation energy balance: SWREB
%--------------------------------------------------------------------------
% SWREBRoofImp		=	Energy Balance shortwave radiation roof impervious area [W/m2 horizontal roof area]
% SWREBRoofVeg		=	Energy Balance shortwave radiation roof vegetated area [W/m2 horizontal roof area]
% SWREBGroundImp	=	Energy Balance shortwave radiation ground impervious area [W/m2 horizontal ground area]
% SWREBGroundBare	=	Energy Balance shortwave radiation ground bare area [W/m2 horizontal ground area]
% SWREBGroundVeg	=	Energy Balance shortwave radiation ground vegetated area [W/m2 horizontal ground area]
% SWREBTree			=	Energy Balance shortwave radiation tree canopy [W/m2 horizontally projected tree area: 4*radius]
% SWREBWallSun		=	Energy Balance shortwave radiation sunlit area [W/m2 vertical wall area]
% SWREBWallShade	=	Energy Balance shortwave radiation shaded area [W/m2 vertical wall area]
% SWREBTotalRoof	=	Energy Balance total shortwave radiation by the roof area [W/m2 horizontal roof area]
% SWREBTotalGround	=	Energy Balance total shortwave radiation by the canyon ground area [W/m2]
% SWREBTotalCanyon	=	Energy Balance total shortwave radiation by all the canyon facets [W/m2]
% SWREBTotalUrban	=	Energy Balance total outgoing shortwave radiation by all the urban elements (roof plus canyon) [W/m2]

SWRabsNames	=	{'SWRabsRoofImp';'SWRabsRoofVeg';'SWRabsTotalRoof';'SWRabsGroundImp';'SWRabsGroundBare';...
					'SWRabsGroundVeg';'SWRabsTree';'SWRabsWallSun';'SWRabsWallShade';...
					'SWRabsTotalGround';'SWRabsTotalCanyon';'SWRabsTotalUrban'};
				
SWRinNames	=	{'SWRinRoofImp';'SWRinRoofVeg';'SWRinTotalRoof';'SWRinGroundImp';'SWRinGroundBare';...
					'SWRinGroundVeg';'SWRinTree';'SWRinWallSun';'SWRinWallShade';...
					'SWRinTotalGround';'SWRinTotalCanyon';'SWRinTotalUrban'};
								
SWRoutNames	=	{'SWRoutRoofImp';'SWRoutRoofVeg';'SWRoutTotalRoof';'SWRoutGroundImp';'SWRoutGroundBare';...
					'SWRoutGroundVeg';'SWRoutTree';'SWRoutWallSun';'SWRoutWallShade';...
					'SWRoutTotalGround';'SWRoutTotalCanyon';'SWRoutTotalUrban'};
				
SWREBNames	=	{'SWREBRoofImp';'SWREBRoofVeg';'SWREBTotalRoof';'SWREBGroundImp';'SWREBGroundBare';...
					'SWREBGroundVeg';'SWREBTree';'SWREBWallSun';'SWREBWallShade';...
					'SWREBTotalGround';'SWREBTotalCanyon';'SWREBTotalUrban'};

for i=1:size(SWRabsNames,1)
	SWRabs.(cell2mat(SWRabsNames(i)))	=	zeros(n,1,m);
	SWRin.(cell2mat(SWRinNames(i)))		=	zeros(n,1,m);
	SWRout.(cell2mat(SWRoutNames(i)))	=	zeros(n,1,m);
	SWREB.(cell2mat(SWREBNames(i)))		=	zeros(n,1,m);
end

%% Longwave radiation
%--------------------------------------------------------------------------
% Absorbed longwave radiation: LWRabs
%--------------------------------------------------------------------------
% LWRabsRoofImp		=	Absorbed longwave radiation roof impervious area [W/m2 horizontal roof area]
% LWRabsRoofVeg		=	Absorbed longwave radiation roof vegetated area [W/m2 horizontal roof area]
% LWRabsGroundImp	=	Absorbed longwave radiation ground impervious area [W/m2 horizontal ground area]
% LWRabsGroundBare	=	Absorbed longwave radiation ground bare area [W/m2 horizontal ground area]
% LWRabsGroundVeg	=	Absorbed longwave radiation ground vegetated area [W/m2 horizontal ground area]
% LWRabsTree		=	Absorbed longwave radiation tree canopy [W/m2 horizontally projected tree area: 4*radius]
% LWRabsWallSun		=	Absorbed longwave radiation sunlit area [W/m2 vertical wall area]
% LWRabsWallShade	=	Absorbed longwave radiation shaded area [W/m2 vertical wall area]
% LWRabsTotalRoof	=	Total absorbed longwave radiation by the roof area [W/m2 horizontal roof area]
% LWRabsTotalGround	=	Total absorbed longwave radiation by the canyon ground area [W/m2]
% LWRabsTotalCanyon	=	Total absorbed longwave radiation by all the canyon facets [W/m2]
% LWRabsTotalUrban	=	Total absorbed longwave radiation by all the urban elements (roof plus canyon) [W/m2]
%--------------------------------------------------------------------------
% Incoming longwave radiation: LWRin
%--------------------------------------------------------------------------
% LWRinRoofImp		=	Incoming longwave radiation roof impervious area [W/m2 horizontal roof area]
% LWRinRoofVeg		=	Incoming longwave radiation roof vegetated area [W/m2 horizontal roof area]
% LWRinGroundImp	=	Incoming longwave radiation ground impervious area [W/m2 horizontal ground area]
% LWRinGroundBare	=	Incoming longwave radiation ground bare area [W/m2 horizontal ground area]
% LWRinGroundVeg	=	Incoming longwave radiation ground vegetated area [W/m2 horizontal ground area]
% LWRinTree			=	Incoming longwave radiation tree canopy [W/m2 horizontally projected tree area: 4*radius]
% LWRinWallSun		=	Incoming longwave radiation sunlit area [W/m2 vertical wall area]
% LWRinWallShade	=	Incoming longwave radiation shaded area [W/m2 vertical wall area]
% LWRinTotalRoof	=	Total incoming longwave radiation by the roof area [W/m2 horizontal roof area]
% LWRinTotalGround	=	Total incoming longwave radiation by the canyon ground area [W/m2]
% LWRinTotalCanyon	=	Total incoming longwave radiation by all the canyon facets [W/m2]
% LWRinTotalUrban	=	Total incoming longwave radiation by all the urban elements (roof plus canyon) [W/m2]
%--------------------------------------------------------------------------
% Outgoing longwave radiation: LWRout
%--------------------------------------------------------------------------
% LWRoutRoofImp		=	Outgoing longwave radiation roof impervious area [W/m2 horizontal roof area]
% LWRoutRoofVeg		=	Outgoing longwave radiation roof vegetated area [W/m2 horizontal roof area]
% LWRoutGroundImp	=	Outgoing longwave radiation ground impervious area [W/m2 horizontal ground area]
% LWRoutGroundBare	=	Outgoing longwave radiation ground bare area [W/m2 horizontal ground area]
% LWRoutGroundVeg	=	Outgoing longwave radiation ground vegetated area [W/m2 horizontal ground area]
% LWRoutTree		=	Outgoing longwave radiation tree canopy [W/m2 horizontally projected tree area: 4*radius]
% LWRoutWallSun		=	Outgoing longwave radiation sunlit area [W/m2 vertical wall area]
% LWRoutWallShade	=	Outgoing longwave radiation shaded area [W/m2 vertical wall area]
% LWRoutTotalRoof	=	Total outgoing longwave radiation by the roof area [W/m2 horizontal roof area]
% LWRoutTotalGround	=	Total outgoing longwave radiation by the canyon ground area [W/m2]
% LWRoutTotalCanyon	=	Total outgoing longwave radiation by all the canyon facets [W/m2]
% LWRoutTotalUrban	=	Total outgoing longwave radiation by all the urban elements (roof plus canyon) [W/m2]
%--------------------------------------------------------------------------
% Energy Balance of longwave radiation: LWREB
%--------------------------------------------------------------------------
% LWREBRoofImp		=	Energy Balance longwave radiation roof impervious area [W/m2 horizontal roof area]
% LWREBRoofVeg		=	Energy Balance longwave radiation roof vegetated area [W/m2 horizontal roof area]
% LWREBGroundImp	=	Energy Balance longwave radiation ground impervious area [W/m2 horizontal ground area]
% LWREBGroundBare	=	Energy Balance longwave radiation ground bare area [W/m2 horizontal ground area]
% LWREBGroundVeg	=	Energy Balance longwave radiation ground vegetated area [W/m2 horizontal ground area]
% LWREBTree			=	Energy Balance longwave radiation tree canopy [W/m2 horizontally projected tree area: 4*radius]
% LWREBWallSun		=	Energy Balance longwave radiation sunlit area [W/m2 vertical wall area]
% LWREBWallShade	=	Energy Balance longwave radiation shaded area [W/m2 vertical wall area]
% LWREBTotalRoof	=	Energy Balance total longwave radiation by the roof area [W/m2 horizontal roof area]
% LWREBTotalGround	=	Energy Balance total longwave radiation by the canyon ground area [W/m2]
% LWREBTotalCanyon	=	Energy Balance total longwave radiation by all the canyon facets [W/m2]
% LWREBTotalUrban	=	Energy Balance total outgoing longwave radiation by all the urban elements (roof plus canyon) [W/m2]

LWRabsNames	=	{'LWRabsRoofImp';'LWRabsRoofVeg';'LWRabsTotalRoof';'LWRabsGroundImp';'LWRabsGroundBare';...
					'LWRabsGroundVeg';'LWRabsTree';'LWRabsWallSun';'LWRabsWallShade';...
					'LWRabsTotalGround';'LWRabsTotalCanyon';'LWRabsTotalUrban'};
				
LWRinNames	=	{'LWRinRoofImp';'LWRinRoofVeg';'LWRinTotalRoof';'LWRinGroundImp';'LWRinGroundBare';...
					'LWRinGroundVeg';'LWRinTree';'LWRinWallSun';'LWRinWallShade';...
					'LWRinTotalGround';'LWRinTotalCanyon';'LWRinTotalUrban'};
								
LWRoutNames	=	{'LWRoutRoofImp';'LWRoutRoofVeg';'LWRoutTotalRoof';'LWRoutGroundImp';'LWRoutGroundBare';...
					'LWRoutGroundVeg';'LWRoutTree';'LWRoutWallSun';'LWRoutWallShade';...
					'LWRoutTotalGround';'LWRoutTotalCanyon';'LWRoutTotalUrban'};
				
LWREBNames	=	{'LWREBRoofImp';'LWREBRoofVeg';'LWREBTotalRoof';'LWREBGroundImp';'LWREBGroundBare';...
					'LWREBGroundVeg';'LWREBTree';'LWREBWallSun';'LWREBWallShade';...
					'LWREBTotalGround';'LWREBTotalCanyon';'LWREBTotalUrban'};

for i=1:size(LWRabsNames,1)
	LWRabs.(cell2mat(LWRabsNames(i)))	=	zeros(n,1,m);
	LWRin.(cell2mat(LWRinNames(i)))		=	zeros(n,1,m);
	LWRout.(cell2mat(LWRoutNames(i)))	=	zeros(n,1,m);
	LWREB.(cell2mat(LWREBNames(i)))		=	zeros(n,1,m);
end

%--------------------------------------------------------------------------
%% Sensible heat flux: Hflux
%--------------------------------------------------------------------------
% HfluxRoofImp		=	Sensible heat flux of impervious roof area to atmosphere [W/m2 horizontal roof area]
% HfluxRoofVeg		=	Sensible heat flux of vegetated roof area to atmosphere [W/m2 horizontal roof area]
% HfluxGroundImp	=	Sensible heat flux of impervious ground area to canyon [W/m2 horizontal ground area]
% HfluxGroundBare	=	Sensible heat flux of bare ground area to canyon [W/m2 horizontal ground area]
% HfluxGroundVeg	=	Sensible heat flux of vegetated ground area to canyon [W/m2 horizontal ground area]
% HfluxGround		=	Sensible heat flux of impground area to canyon [W/m2 horizontal ground area]
% HfluxTree			=	Sensible heat flux of tree canopy to canyon [W/m2 horizontally projected tree area: 4*radius]
% HfluxWallSun		=	Sensible heat flux of sunlit wall to canyon [W/m2 vertical wall area]
% HfluxWallShade	=	Sensible heat flux of shaded wall to canyon [W/m2 vertical wall area]
% HfluxCanyon		=	Sensible heat flux of canyon to atmosphere [W/m2 horizontal canyon area]
% HfluxRoof			=	Total sensible heat flux of roof area to atmosphere [W/m2 horizontal roof area]
% HfluxUrban		=	Total sensible heat flux of urban area to atmosphere [W/m2 horizontal area]
				
HfluxNames	=	{'HfluxRoofImp';'HfluxRoofVeg';'HfluxRoof';'HfluxGroundImp';'HfluxGroundBare';...
					'HfluxGroundVeg';'HfluxGround';'HfluxTree';'HfluxWallSun';'HfluxWallShade';...
					'HfluxCanyon';'HfluxUrban'};

for i=1:size(HfluxNames,1)
	Hflux.(cell2mat(HfluxNames(i)))	=	zeros(n,1,m);
end

%--------------------------------------------------------------------------
%% Latent heat flux: LEflux
%--------------------------------------------------------------------------
% LEfluxRoofImp			=	Latent heat flux of intercepted water from impervious roof area to atmosphere (LEroof_imp_pond) [W/m2 horizontal roof area]
% LEfluxRoofVegInt		=	Latent heat flux of intercepted water on roof vegetation to atmosphere (LEroof_veg_int) [W/m2 horizontal roof area]
% LEfluxRoofVegPond		=	Latent heat flux of intercepted water on ground under roof vegetation to atmosphere (LEroof_veg_pond) [W/m2 horizontal roof area]
% LEfluxRoofVegSoil		=	Latent heat flux of water from roof soil under vegetation to atmosphere (LEroof_veg_soil) [W/m2 horizontal roof area]
% LTEfluxRoofVeg		=	Latent heat flux of transpiration from roof plants to atmosphere (LTEroof_veg) [W/m2 horizontal roof area]
% LEfluxRoofVeg			=	Total latent heat flux of vegetated roof to atmosphere [W/m2 horizontal roof area]
% LEfluxRoof			=	Total latent heat flux of roof to atmosphere [W/m2 horizontal roof area]
% LEfluxGroundImp		=	Latent heat flux of intercepted water on impervious ground area to canyon (LEground_imp_pond)[W/m2 horizontal ground area]
% LEfluxGroundBarePond	=	Latent heat flux of  water on bare ground to canyon (LEground_bare_pond)[W/m2 horizontal ground area]
% LEfluxGroundBareSoil	=	Latent heat flux of  water from bare ground to canyon (LEground_bare_soil) [W/m2 horizontal ground area]
% LEfluxGroundBare		=	Total latent heat flux of bare ground area to canyon [W/m2 horizontal ground area]
% LEfluxGroundVegInt	=	Latent heat flux of intercepted water on ground vegetation to canyon (LEground_veg_int) [W/m2 horizontal ground area]
% LEfluxGroundVegPond	=	Latent heat flux of intercepted water on ground under vegetation to canyon (LEground_veg_pond) [W/m2 horizontal ground area]
% LEfluxGroundVegSoil	=	Latent heat flux of water from ground soil under vegetation to canyon (LEground_veg_soil) [W/m2 horizontal ground area]
% LTEfluxGroundVeg		=	Latent heat flux of transpiration from ground plants to canyon (LTEground_veg) [W/m2 horizontal ground area]
% LEfluxGroundVeg		=	Total latent heat flux of vegetated ground to canyon [W/m2 horizontal ground area]
% LEfluxGround			=	Total latent heat flux of ground to canyon [W/m2 horizontal roof area]
% LEfluxTreeInt			=	Latent heat flux of intercepted water on tree canopy to canyon (LE_tree_int) [W/m2 horizontally projected tree area: 4*radius]
% LTEfluxTree			=	Latent heat flux of transpiration from tree canopy to canyon (LTE_tree) [W/m2 horizontally projected tree area: 4*radius]
% LEfluxTree			=	Total latent heat flux of tree canopy to canyon [W/m2 horizontally projected tree area: 4*radius]
% LEfluxWallSun			=	Latent heat flux of sunlit wall to canyon [W/m2 vertical wall area]
% LEfluxWallShade		=	Latent heat flux of shaded wall to canyon [W/m2 vertical wall area]
% LEfluxCanyon			=	Latent heat flux of canyon to atmosphere [W/m2 horizontal canyon area]
% LEfluxUrban			=	Total latent heat flux of urban area to atmosphere [W/m2 horizontal area]

LEfluxNames	=	{'LEfluxRoofImp';'LEfluxRoofVegInt';'LEfluxRoofVegPond';'LEfluxRoofVegSoil';...
					'LTEfluxRoofVeg';'LEfluxRoofVeg';'LEfluxRoof';'LEfluxGroundImp';'LEfluxGroundBarePond';...
					'LEfluxGroundBareSoil';'LEfluxGroundBare';'LEfluxGroundVegInt';...
					'LEfluxGroundVegPond';'LEfluxGroundVegSoil';'LTEfluxGroundVeg';'LEfluxGroundVeg';...
					'LEfluxGround';'LEfluxTreeInt';'LTEfluxTree';'LEfluxTree';'LEfluxWallSun';...
					'LEfluxWallShade';'LEfluxCanyon';'LEfluxUrban'};

for i=1:size(LEfluxNames,1)
	LEflux.(cell2mat(LEfluxNames(i)))	=	zeros(n,1,m);
end

%--------------------------------------------------------------------------
%% Conductive heat fluxes: Gflux
%--------------------------------------------------------------------------
% G1RoofImp		=	Conductive heat flux of first layer of impervious roof [W/m2]
% G1RoofVeg		=	Conductive heat flux of first layer of vegetated roof [W/m2]
% G2RoofImp		=	Conductive heat flux of second layer of impervious roof [W/m2]
% G2RoofVeg		=	Conductive heat flux of second layer of vegetated roof [W/m2]
% G1GroundImp	=	Conductive heat flux of impervious ground (G_groundimp) [W/m2]
% G1GroundBare	=	Conductive heat flux of bare ground (G_groundbare) [W/m2]
% G1GroundVeg	=	Conductive heat flux of vegetated ground (G_groundveg) [W/m2]
% GTree			=	Conductive heat flux tree [W/m2]
% G1WallSun		=	Conductive heat flux of first layer of sunlit wall (G1_wallsun) [W/m2]
% G1WallShade	=	Conductive heat flux of first layer of shaded wall (G1_wallshade) [W/m2]
% G2WallSun		=	Conductive heat flux of second layer of sunlit wall (G2_wallsun) [W/m2]
% G2WallShade	=	Conductive heat flux of second layer of shaded wall (G2_wallshade) [W/m2]
	
GfluxNames	=	{'G1RoofImp';'G1RoofVeg';'G2RoofImp';'G2RoofVeg';'G1Roof';'G2Roof';...
					'G1GroundImp';'G1GroundBare';'G1GroundVeg';'G1Ground';'GTree';'G1WallSun';...
					'G1WallShade';'G2WallSun';'G2WallShade';'G1Canyon';'G2Canyon';'G1Urban';'G2Urban'};

for i=1:size(GfluxNames,1)
	Gflux.(cell2mat(GfluxNames(i)))	=	zeros(n,1,m);
end

%--------------------------------------------------------------------------
%% Heat storage in surfaces: dStorage
%--------------------------------------------------------------------------
% dsRoofImp		=	Storage of energy in impervious roof [W/m2]
% dsRoofVeg		=	Storage of energy in vegetated roof [W/m2]
% dsGroundImp	=	Storage of energy in impervious ground [W/m2]
% dsGroundBare	=	Storage of energy in bare ground [W/m2]
% dsGroundVeg	=	Storage of energy in vegetated ground [W/m2]
% dsTree		=	Storage of energy in tree canopy [W/m2]
% dsWallSun		=	Storage of energy in sunlit wall  [W/m2]
% dsWallShade	=	Storage of energy in shaded wall [W/m2]
% dsCanyonAir	=	Storage of energy in canyon air [W/m2]

dStorageNames	=	{'dsRoofImp';'dsRoofVeg';'dsRoof';'dsGroundImp';'dsGroundBare';...
					'dsGroundVeg';'dsTree';'dsWallSun';'dsWallShade';'dsCanyonAir'};

for i=1:size(dStorageNames,1)
	dStorage.(cell2mat(dStorageNames(i)))	=	zeros(n,1,m);
end

%--------------------------------------------------------------------------
%% Resistances: RES
%--------------------------------------------------------------------------
% raRooftoAtm	=	Aerodynamic resistance ra from roof to atmosphere [s/m]
% rap_LRoof		=	Undercanopy resistance rap_L roof [s/m]
% rb_LRoof		=	Leaf boundary resistance rb_L roof [s/m]
% r_soilRoof	=	Soil resistance rb_soil roof [s/m]
% rs_sunRoof	=	Stomata resistance sunlit vegetation rs_sun_roof [s/m]
% rs_shdRoof	=	Stomata resistance shaded vegetation rs_shd_roof [s/m]
% raCanyontoAtm	=	Aerodynamic resistance ra from canyon to atmosphere [s/m]
% raGroundtoAtm	=	Aerodynamic resistance ra from ground to canyon air [s/m]
% raTreetoAtm	=	Aerodynamic resistance ra from tree to canyon air [s/m]
% raWalltoAtm	=	Aerodynamic resistance ra from wall to canyon air [s/m]
% rap_HGround	=	Undercanopy resistance rap_H ground [s/m]
% rap_LGround	=	Undercanopy resistance rap_L ground [s/m]
% rb_HGround	=	Leaf boundary resistance rb_H ground [s/m]
% rb_LGround	=	Leaf boundary resistance rb_L ground [s/m]
% r_soilGround	=	Soil resistance rb_soil ground [s/m]
% rs_sunGround	=	Stomata resistance sunlit vegetation rs_sun_ground [s/m]
% rs_shdGround	=	Stomata resistance shaded vegetation rs_shd_ground [s/m]
% rs_sunTree	=	Stomata resistance sunlit vegetation rs_sun_tree [s/m]
% rs_shdTree	=	Stomata resistance shaded vegetation rs_shd_ground [s/m]
	
RESNames	=	{'raRooftoAtm';'rap_LRoof';'rb_LRoof';'r_soilRoof';'rs_sunRoof';'rs_shdRoof';...
					'raCanyontoAtm';'rap_can';'rap_Htree_In';'rb_HGround';'rb_LGround';...
					'r_soilGroundbare';'r_soilGroundveg';'alp_soilGroundbare';'alp_soilGroundveg';...
					'rs_sunGround';'rs_shdGround';'rs_sunTree';'rs_shdTree';...
					'RES_w1';'RES_w2';'rap_W1_In';'rap_W2_In';'rap_Zp1'};

for i=1:size(RESNames,1)
	RES.(cell2mat(RESNames(i)))	=	zeros(n,1,m);
end

%--------------------------------------------------------------------------
%% Water fluxes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% Evapotranspiration: Eflux
%--------------------------------------------------------------------------
% EfluxRoofImp			=	Evaporation flux of intercepted water from impervious roof area to atmosphere (Eroof_imp_pond) [kg/m^2*s horizontal roof area]
% EfluxRoofVegInt		=	Evaporation flux of intercepted water on roof vegetation to atmosphere (Eroof_veg_int) [kg/m^2*s horizontal roof area]
% EfluxRoofVegPond		=	Evaporation flux of intercepted water on ground under roof vegetation to atmosphere (Eroof_veg_pond) [kg/m^2*s horizontal roof area]
% EfluxRoofVegSoil		=	Evaporation flux of water from roof soil under vegetation to atmosphere (Eroof_veg_soil) [kg/m^2*s horizontal roof area]
% TEfluxRoofVeg			=	Evaporation flux of transpiration from roof plants to atmosphere (TEroof_veg) [kg/m^2*s horizontal roof area]
% EfluxRoofVeg			=	Total evaporation flux of vegetated roof to atmosphere [kg/m^2*s horizontal roof area]
% EfluxRoof				=	Total evaporation flux of roof to atmosphere [kg/m^2*s horizontal roof area]
% EfluxGroundImp		=	Evaporation flux of intercepted water on impervious ground area to canyon (Eground_imp_pond)[kg/m^2*s horizontal ground area]
% EfluxGroundBarePond	=	Evaporation flux of  water on bare ground to canyon (Eground_bare_pond)[kg/m^2*s horizontal ground area]
% EfluxGroundBareSoil	=	Evaporation flux of  water from bare ground to canyon (Eground_bare_soil) [kg/m^2*s horizontal ground area]
% EfluxGroundBare		=	Total evaporation flux of bare ground area to canyon [kg/m^2*s horizontal ground area]
% EfluxGroundVegInt		=	Evaporation flux of intercepted water on ground vegetation to canyon (Eground_veg_int) [kg/m^2*s horizontal ground area]
% EfluxGroundVegPond	=	Evaporation flux of intercepted water on ground under vegetation to canyon (Eground_veg_pond) [kg/m^2*s horizontal ground area]
% EfluxGroundVegSoil	=	Evaporation flux of water from ground soil under vegetation to canyon (Eground_veg_soil) [kg/m^2*s horizontal ground area]
% TEfluxGroundVeg		=	Evaporation flux of transpiration from ground plants to canyon (TEground_veg) [kg/m^2*s horizontal ground area]
% EfluxGroundVeg		=	Total evaporation flux of vegetated ground to canyon [kg/m^2*s horizontal ground area]
% EfluxGround			=	Total evaporation flux of ground to canyon [kg/m^2*s horizontal roof area]
% EfluxTreeInt			=	Evaporation flux of intercepted water on tree canopy to canyon (E_tree_int) [kg/m^2*s horizontally projected tree area: 4*radius]
% TEfluxTree			=	Evaporation flux of transpiration from tree canopy to canyon (TE_tree) [kg/m^2*s horizontally projected tree area: 4*radius]
% EfluxTree				=	Total evaporation flux of tree canopy to canyon [kg/m^2*s horizontally projected tree area: 4*radius]
% EfluxWallSun			=	Evaporation flux of sunlit wall to canyon [kg/m^2*s vertical wall area]
% EfluxWallShade		=	Evaporation flux of shaded wall to canyon [kg/m^2*s vertical wall area]
% EfluxCanyon			=	Evaporation flux of canyon to atmosphere [kg/m^2*s horizontal canyon area]
% EfluxUrban			=	Total evaporation flux of urban area to atmosphere kg/m^2*s horizontal area]

EfluxNames	=	{'EfluxRoofImp';'EfluxRoofVegInt';'EfluxRoofVegPond';'EfluxRoofVegSoil';...
					'TEfluxRoofVeg';'EfluxRoofVeg';'EfluxRoof';'EfluxGroundImp';'EfluxGroundBarePond';...
					'EfluxGroundBareSoil';'EfluxGroundBare';'EfluxGroundVegInt';...
					'EfluxGroundVegPond';'EfluxGroundVegSoil';'TEfluxGroundVeg';'EfluxGroundVeg';...
					'EfluxGround';'EfluxTreeInt';'TEfluxTree';'EfluxTree';'EfluxWallSun';...
					'EfluxWallShade';'EfluxCanyon';'EfluxUrban'};

for i=1:size(EfluxNames,1)
	Eflux.(cell2mat(EfluxNames(i)))	=	zeros(n,1,m);
end

%--------------------------------------------------------------------------
%% Runoff: Runoff and Runon: Runon
%--------------------------------------------------------------------------
% QRoofImp			=	Runoff of impervious area of roof (q_runon_imp) [mm/time step] 
% QRoofVegDrip		=	Runoff, Dripping, etc from ground vegetation to roof ground (q_runon_veg) [mm/time step] 
% QRoofVegPond		=	Runoff from roof ground under vegetation due to limitation in infiltration capacity (q_runon_ground_veg) [mm/time step] 
% QRoofVegSoil		=	Runoff due to roof soil saturation (Rd_veg)[mm/time step] 
% QGroundImp		=	Runoff of impervious area of ground (q_runon_imp)[mm/time step] 
% QGroundBarePond	=	Runoff of bare area of ground due to limitation in infiltration capacity(q_runon_bare)[mm/time step] 
% QGroundBareSoil	=	Runoff of bare area of ground due to soil saturation (Rd_bare)[mm/time step] 
% QTree				=	Runoff, Dripping, etc from tree to ground (q_runon_tree)[mm/time step] 
% QGroundVegDrip	=	Runoff, Dripping, etc from ground vegetation to roof ground (q_runon_veg)[mm/time step] 
% QGroundVegPond	=	Runoff from roof ground under vegetation due to limitation in infiltration capacity (q_runon_ground_veg)[mm/time step] 
% QGroundVegSoil	=	Runoff due to roof soil saturation (Rd_veg)

RunoffNames	=	{'QRoofImp';'QRoofVegDrip';'QRoofVegPond';'QRoofVegSoil';...
					'QGroundImp';'QGroundBarePond';'QGroundBareSoil';'QTree';'QGroundVegDrip';...
					'QGroundVegPond';'QGroundVegSoil'};

for i=1:size(RunoffNames,1)
	Runoff.(cell2mat(RunoffNames(i)))	=	zeros(n,1,m);
end

% RunonRoofTot		=	Total roof runon to the next time step [mm/time step] 
% RunoffRoofTot		=	Total roof runoff that is removed from the system [mm/time step] 
% RunonGroundTot	=	Total ground runon to the next time step [mm/time step] 
% RunoffGroundTot	=	Total ground runoff that is removed from the system [mm/time step] 

RunonNames	=	{'RunonRoofTot';'RunoffRoofTot';'RunonGroundTot';'RunoffGroundTot';'RunonUrban';'RunoffUrban'};

for i=1:size(RunonNames,1)
	Runon.(cell2mat(RunonNames(i)))	=	zeros(n,1,m);
end

%--------------------------------------------------------------------------
%% Leakage: Leakage
%--------------------------------------------------------------------------
% LkRoofImp		=	Leakage from impervious roof (Lk_imp)[mm/h]
% LkRoofVeg		=	Leakage from last soil layer of vegetated roof (Lk_soil_veg)[mm/h]
% LkRoof		=	Total leakage of roof [mm/h]
% LkGroundImp	=	Leakage from impervious ground (Lk_imp)[mm/h]
% LkGroundBare	=	Leakage from last soil layer of bare ground (Lk_soil_bare)[mm/h]
% LkGroundVeg	=	Leakage from last soil layer of vegetated ground (Lk_soil_veg)[mm/h]
% LkGround		=	Total leakage of ground[mm/h]

LeakageNames	=	{'LkRoofImp';'LkRoofVeg';'LkRoof';'LkGroundImp';...
					'LkGroundBare';'LkGroundVeg';'LkGround';'LkUrban'};

for i=1:size(LeakageNames,1)
	Leakage.(cell2mat(LeakageNames(i)))	=	zeros(n,1,m);
end

%--------------------------------------------------------------------------
%% Interception: Int
%--------------------------------------------------------------------------
% IntRoofImp		=	Interception on impervious roof area(In_ground_imp) [mm]
% IntRoofVegPlant	=	Interception on plant surfaces (In_ground_veg) [mm]
% IntRoofVegGround	=	Interception on ground (In_ground_underveg) [mm]
% IntGroundImp		=	Interception on impervious ground area (In_ground_imp) [mm]
% IntGroundBare		=	Interception on bare ground area (In_ground_bare) [mm]
% IntGroundVegPlant	=	Interception on plant surfaces (In_ground_veg) [mm]
% IntGroundVegGround=	Interception on ground (In_ground_underveg) [mm]
% IntTree			=	Interception on tree (In_tree) [mm]

IntNames	=	{'IntRoofImp';'IntRoofVegPlant';'IntRoofVegGround';'IntRooftot';'IntGroundImp';...
					'IntGroundBare';'IntGroundVegPlant';'IntGroundVegGround';'IntTree'};

for i=1:size(IntNames,1)
	Int.(cell2mat(IntNames(i)))	=	zeros(n,1,m);
end

%--------------------------------------------------------------------------
% Change in interception: dInt_dt
%--------------------------------------------------------------------------
% dInt_dtRoofImp		=	Change in interception on impervious roof area (dIn_imp_dt)[mm/h]
% dInt_dtRoofVegPlant	=	Change in interception on plant surfaces (dIn_veg_dt)[mm/h]
% dInt_dtRoofVegGround	=	Change in interception on ground (dIn_ground_veg_dt)[mm/h]
% dInt_dtGroundImp		=	Change in interception on impervious ground area (dIn_imp_dt)[mm/h]
% dInt_dtGroundBare		=	Change in interception on bare ground area (dIn_bare_dt)[mm/h]
% dInt_dtGroundVegPlant	=	Change in interception on plant surfaces (dIn_veg_dt)[mm/h]
% dInt_dtGroundVegGround=	Change in interception on ground (dIn_ground_veg_dt)[mm/h]
% dInt_dtTree			=	Change in interception on tree (dIn_tree_dt)[mm/h]

dInt_dtNames	=	{'dInt_dtRoofImp';'dInt_dtRoofVegPlant';'dInt_dtRoofVegGround';'dInt_dtRooftot';'dInt_dtGroundImp';...
					'dInt_dtGroundBare';'dInt_dtGroundVegPlant';'dInt_dtGroundVegGround';'dInt_dtTree'};

for i=1:size(dInt_dtNames,1)
	dInt_dt.(cell2mat(dInt_dtNames(i)))	=	zeros(n,1,m);
end

%--------------------------------------------------------------------------
%% Infiltration: Infiltration
%--------------------------------------------------------------------------
% fRoofVeg		=	Infiltration in first soil layer of vegetated roof (f_roof_veg)[mm/h]
% fGroundBare	=	Infiltration in first soil layer of bare ground (f_ground_bare)[mm/h]
% fGroundVeg	=	Infiltration in first soil layer of vegetated ground (f_ground_veg)	[mm/h]

InfiltrationNames	=	{'fRoofVeg';'fGroundBare';'fGroundVeg';'fGroundImp'};

for i=1:size(InfiltrationNames,1)
	Infiltration.(cell2mat(InfiltrationNames(i)))	=	zeros(n,1,m);
end

%--------------------------------------------------------------------------
%% Water volume in soil: dVwater_dt and Change in water volume in soil: dVwater_dt
%--------------------------------------------------------------------------
% Initializing soil water content in the first time step.
% I chose field capacity O33 as a starting point
Vwater						=	[];
Vwater.VRoofSoilVeg			=	zeros(n,ParSoil.Roof.ms,m);		%  Water volume in the different soil layers of roof (Vw_soil) [mm]
Vwater.VGroundSoilImp		=	zeros(n,ParSoil.Ground.ms,m);	%  Water volume in the different soil layers of ground (Vw_soil)[mm]
Vwater.VGroundSoilBare		=	zeros(n,ParSoil.Ground.ms,m);	%  Water volume in the different soil layers of ground (Vw_soil)[mm]
Vwater.VGroundSoilVeg		=	zeros(n,ParSoil.Ground.ms,m);	%  Water volume in the different soil layers of ground (Vw_soil)[mm]
Vwater.VGroundSoilTot		=	zeros(n,ParSoil.Ground.ms,m);	%  Water volume in the different soil layers of ground (Vw_soil)[mm]

Vwater.VRoofSoilVeg(1,:,:)	=	repmat(ParSoil.Roof.O33.*ParSoil.Roof.dz,1,1,m);		% Starting point at field capacity[mm]
Vwater.VGroundSoilImp(1,:,:)=	repmat(ParSoil.Ground.O33.*ParSoil.Ground.dz,1,1,m);	% Starting point at field capacity[mm]
Vwater.VGroundSoilBare(1,:,:)=	repmat(ParSoil.Ground.O33.*ParSoil.Ground.dz,1,1,m);	% Starting point at field capacity[mm]
Vwater.VGroundSoilVeg(1,:,:)=	repmat(ParSoil.Ground.O33.*ParSoil.Ground.dz,1,1,m);	% Starting point at field capacity[mm]
Vwater.VGroundSoilTot(1,:,:)=	repmat(ParSoil.Ground.O33.*ParSoil.Ground.dz,1,1,m);	% Starting point at field capacity[mm]

dVwater_dt						=	[];
dVwater_dt.dVRoofSoilVeg_dt		=	zeros(n,1,m);	% Maybe only total? I should check zeros(n,1);[mm]
dVwater_dt.dVGroundSoilImp_dt	=	zeros(n,1,m); % Maybe only total? I should check zeros(n,1);[mm]
dVwater_dt.dVGroundSoilBare_dt	=	zeros(n,1,m); % Maybe only total? I should check zeros(n,1);[mm]
dVwater_dt.dVGroundSoilVeg_dt	=	zeros(n,1,m); % Maybe only total? I should check zeros(n,1);[mm]
dVwater_dt.dVGroundSoilTot_dt	=	zeros(n,1,m); % Maybe only total? I should check zeros(n,1);[mm]

%--------------------------------------------------------------------------
%% Soil moisture: Owater
%--------------------------------------------------------------------------
Owater						=	[];
Owater.OwRoofSoilVeg		=	zeros(n,ParSoil.Roof.ms,m);			%  Soil moisture in the different soil layers of roof [-]
Owater.OwGroundSoilImp		=	zeros(n,ParSoil.Ground.ms,m);		%  Soil moisture in the different soil layers of ground [-]
Owater.OwGroundSoilBare		=	zeros(n,ParSoil.Ground.ms,m);		%  Soil moisture in the different soil layers of ground [-]
Owater.OwGroundSoilVeg		=	zeros(n,ParSoil.Ground.ms,m);		%  Soil moisture in the different soil layers of ground [-]
Owater.OwGroundSoilTot		=	zeros(n,ParSoil.Ground.ms,m);		%  Soil moisture in the different soil layers of ground [-]

Owater.OwRoofSoilVeg(1,:,:)	=	ParSoil.Roof.O33;				% Starting point at field capacity
Owater.OwGroundSoilImp(1,:,:)=	ParSoil.Ground.O33;				% Starting point at field capacity
Owater.OwGroundSoilBare(1,:,:)=	ParSoil.Ground.O33;				% Starting point at field capacity
Owater.OwGroundSoilVeg(1,:,:)=	ParSoil.Ground.O33;				% Starting point at field capacity
Owater.OwGroundSoilTot(1,:,:)=	ParSoil.Ground.O33;				% Starting point at field capacity

Owater.OwGroundSoilImp(:,1:2,:)=	NaN;				% Starting point at field capacity

OSwater						=	[];
OSwater.OSwRoofSoilVeg		=	zeros(n,ParSoil.Roof.ms,m);
OSwater.OSwGroundSoilImp	=	zeros(n,ParSoil.Ground.ms,m);
OSwater.OSwGroundSoilBare	=	zeros(n,ParSoil.Ground.ms,m);
OSwater.OSwGroundSoilVeg	=	zeros(n,ParSoil.Ground.ms,m);
OSwater.OSwGroundSoilTot	=	zeros(n,ParSoil.Ground.ms,m);

%--------------------------------------------------------------------------
% Lateral soil water flux
%--------------------------------------------------------------------------
Qinlat				=	[];
Qinlat.Qin_bare2imp	=	zeros(n,ParSoil.Ground.ms,m);
Qinlat.Qin_veg2imp	=	zeros(n,ParSoil.Ground.ms,m);
Qinlat.Qin_veg2bare	=	zeros(n,ParSoil.Ground.ms,m);
Qinlat.Qin_imp2bare	=	zeros(n,ParSoil.Ground.ms,m);
Qinlat.Qin_bare2veg	=	zeros(n,ParSoil.Ground.ms,m);
Qinlat.Qin_imp2veg	=	zeros(n,ParSoil.Ground.ms,m);
Qinlat.Qin_imp		=	zeros(n,ParSoil.Ground.ms,m);
Qinlat.Qin_bare		=	zeros(n,ParSoil.Ground.ms,m);
Qinlat.Qin_veg		=	zeros(n,ParSoil.Ground.ms,m);

%% Max extractable water and soil water potential for plants in soil
%--------------------------------------------------------------------------
% Max extractable water: ExWater
%--------------------------------------------------------------------------
% Extractable water: Soil moisture in the different soil layers (Exwat) [mm m2 / m2 ground h ]
ExWaterNames	=	{'ExWaterRoofVeg_H';'ExWaterRoofVeg_L';...
					'ExWaterGroundImp_H';'ExWaterGroundImp_L';...
					'ExWaterGroundBare_H';'ExWaterGroundBare_L';...
					'ExWaterGroundVeg_H';'ExWaterGroundVeg_L';...
					'ExWaterGroundTot_H';'ExWaterGroundTot_L'};
for i=1:2
	ExWater.(cell2mat(ExWaterNames(i)))	=	zeros(n,ParSoil.Roof.ms,m);
end

for i=3:size(ExWaterNames,1)
	ExWater.(cell2mat(ExWaterNames(i)))	=	zeros(n,ParSoil.Ground.ms,m);
end

%--------------------------------------------------------------------------
% Soil water potential for plants in soil: SoilPotW
%--------------------------------------------------------------------------
% SoilPotWRoof_H	=	soil water potential for plants (Psi_s_H)[MPa]
% SoilPotWRoof_L	=	soil water potential for plants (Psi_s_L)[MPa]
% SoilPotWGround_H	=	soil water potential for plants (Psi_s_H)[MPa]
% SoilPotWGround_L	=	soil water potential for plants (Psi_s_L)[MPa]
SoilPotWNames	=	{'SoilPotWRoofVeg_H';'SoilPotWRoofVeg_L';...
					'SoilPotWGroundImp_H';'SoilPotWGroundImp_L';...
					'SoilPotWGroundBare_H';'SoilPotWGroundBare_L';...
					'SoilPotWGroundVeg_H';'SoilPotWGroundVeg_L';...
					'SoilPotWGroundTot_H';'SoilPotWGroundTot_L'};

for i=1:size(SoilPotWNames,1)
	SoilPotW.(cell2mat(SoilPotWNames(i)))	=	zeros(n,1,m);
end

%--------------------------------------------------------------------------
%% Ci = Leaf Interior  CO2 concentration: CiCO2Leaf 
%--------------------------------------------------------------------------
% CiCO2LeafGroundVegSun	=	Ci_sun_veg [umolCO2/mol]
% CiCO2LeafGroundVegShd	=	Ci_shd_veg [umolCO2/mol]
% CiCO2LeafTreeSun		=	Ci_sun_tree [umolCO2/mol]
% CiCO2LeafTreeShd		=	Ci_shd_tree [umolCO2/mol]

CiCO2LeafNames	=	{'CiCO2LeafRoofVegSun';'CiCO2LeafRoofVegShd';...
					'CiCO2LeafGroundVegSun';'CiCO2LeafGroundVegShd';...
					'CiCO2LeafTreeSun';'CiCO2LeafTreeShd'};

for i=1:size(CiCO2LeafNames,1)
	CiCO2Leaf.(cell2mat(CiCO2LeafNames(i)))			=	zeros(n,1,m);
	CiCO2Leaf.(cell2mat(CiCO2LeafNames(i)))(1,:,:)	=	400;
end

%--------------------------------------------------------------------------
%% Energy and Water Balance
%--------------------------------------------------------------------------
WBRoofNames	=	{'WBRoofImp';'WBRoofVegInVeg';'WBRoofVegInGround';'WBRoofVegSoil';...
				'WBRoofVeg';'WBRoofTot'};

for i=1:size(WBRoofNames,1)
	WBRoof.(cell2mat(WBRoofNames(i)))	=	zeros(n,1,m);
end


WBCanyonIndvNames	=	{'WB_In_tree';'WB_In_gveg';'WB_In_gimp';'WB_In_gbare';...
						'WB_Pond_gveg';'WB_Soil_gimp';'WB_Soil_gbare';'WB_Soil_gveg'};

for i=1:size(WBCanyonIndvNames,1)
	WBCanyonIndv.(cell2mat(WBCanyonIndvNames(i)))	=	zeros(n,1,m);
end


WBCanyonTotNames	=	{'WBsurf_tree';'WBsurf_imp';'WBsurf_bare';'WBsurf_veg';...
						'WBsoil_imp';'WBsoil_bare';'WBsoil_veg';'WBimp_tot';'WBbare_tot';'WBveg_tot';...
						'WBcanyon_flux';'WBtree_level';'WBground_level';'WBsoil_level';'WBcanyon_level'};

for i=1:size(WBCanyonTotNames,1)
	WBCanyonTot.(cell2mat(WBCanyonTotNames(i)))	=	zeros(n,1,m);
end

%--------------------------------------------------------------------------
% Energy Balance
%--------------------------------------------------------------------------
% EBGroundImp	=	EBalance_groundimp [W/m^2]
% EBGroundBare	=	EBalance_groundbare [W/m^2]
% EBGroundVeg	=	EBalance_groundveg [W/m^2]
% EBTree		=	EBalance_tree [W/m^2]
% EBWallSun		=	EBalance_wallsun [W/m^2]
% EBWallShade	=	EBalance_wallshade [W/m^2]
% EBWallSunInt	=	EBalance_wallsun_interior [W/m^2]
% EBWallShadeInt=	EBalance_wallshade_interior [W/m^2]
% EBCanyonT		=	EBalance_canyon_temp [W/m^2]
% EBCanyonQ		=	EBalance_canyon_humid [kg/kg]

EBNames	=	{'EBRoofImp';'EBRoofVeg';'EBGroundImp';'EBGroundBare';...
					'EBGroundVeg';'EBTree';'EBWallSun';'EBWallShade';'EBWallSunInt';...
					'EBWallShadeInt';'EBCanyonT';'EBCanyonQ'};

for i=1:size(EBNames,1)
	EB.(cell2mat(EBNames(i)))	=	zeros(n,1,m);
end

%--------------------------------------------------------------------------
%% Wind speed: Wind
%--------------------------------------------------------------------------
% u_Hcan: Wind speed at canyon calculation height (hdisp + canyon roughness height) (m/s)
% u_Zref_und: Wind speed at undercanopy reference height (m/s)
% u_ZPerson: Wind speed at person height (or point which was specified by the user (m/s)

WindNames	=	{'u_Hcan';'u_Zref_und';'u_ZPerson'};

for i=1:size(WindNames,1)
	Wind.(cell2mat(WindNames(i)))	=	zeros(n,1,m);
end

%--------------------------------------------------------------------------
%% Success of energy balance solver: Solver
%--------------------------------------------------------------------------
% SuccessRoof: Bolean indicateing convergence of solution for roof
% SuccessCanyon: Bolean indicateing convergence of solution for canyon
% ValuesRoof: Energy balance closure for the different roof facets (W/m^2)
% ValuesCanyon: Energy balance closure for the different canyon facets (W/m^2)
% TsolverRoof: Temperatures of different roof facets (K)
% TsolverCanyon: Temperatures and humidity of different canyon factes and air (K), (kg/kg)


Solver				=	[];
Solver.SuccessRoof	=	zeros(n,1,m);
Solver.SuccessCanyon=	zeros(n,1,m);
Solver.ValuesRoof	=	zeros(n,4,m);
Solver.ValuesCanyon	=	zeros(n,10,m);
Solver.TsolverRoof	=	zeros(n,4,m);
Solver.TsolverCanyon=	zeros(n,10,m);

%--------------------------------------------------------------------------
%% Temperature and humidity at 2m canyon height: Results2m
%--------------------------------------------------------------------------
% T2m: 2m air temperature (K)
% q2m: 2m specific humidity (kg/kg)
% e_T2m: 2m vapor pressure (Pa)
% RH_T2m: 2m relative humidity (-)

Results2m		=	[];
Results2m.T2m	=	zeros(n,1,m);
Results2m.q2m	=	zeros(n,1,m);
Results2m.e_T2m	=	zeros(n,1,m);
Results2m.RH_T2m=	zeros(n,1,m);
Results2m.qcan	=	zeros(n,1,m);
Results2m.e_Tcan=	zeros(n,1,m);
Results2m.RH_Tcan=	zeros(n,1,m);

%--------------------------------------------------------------------------
%% Energy fluxes at 2m canyon height: Results2mEnergyFlux
%--------------------------------------------------------------------------
Results2mEnergyFlux				=	[];
Results2mEnergyFlux.DHi			=	zeros(n,1,m);
Results2mEnergyFlux.Himp_2m		=	zeros(n,1,m);
Results2mEnergyFlux.Hbare_2m	=	zeros(n,1,m);
Results2mEnergyFlux.Hveg_2m		=	zeros(n,1,m);
Results2mEnergyFlux.Hwsun_2m	=	zeros(n,1,m);
Results2mEnergyFlux.Hwshade_2m	=	zeros(n,1,m);
Results2mEnergyFlux.Hcan_2m		=	zeros(n,1,m);

Results2mEnergyFlux.DEi			=	zeros(n,1,m);
Results2mEnergyFlux.Eimp_2m		=	zeros(n,1,m);
Results2mEnergyFlux.Ebare_soil_2m=	zeros(n,1,m);
Results2mEnergyFlux.Eveg_int_2m	=	zeros(n,1,m);
Results2mEnergyFlux.Eveg_soil_2m=	zeros(n,1,m);
Results2mEnergyFlux.TEveg_2m	=	zeros(n,1,m);
Results2mEnergyFlux.Ecan_2m		=	zeros(n,1,m);

%--------------------------------------------------------------------------
%% Mean radiant temperature variables: MeanRadiantTemperature
%--------------------------------------------------------------------------
% Tmrt: Mean radiant temperature (deg C)
% BoleanInSun: Point of Tmrt calculation is in sun or in shade
% SWRdir_Person: Direct shortwave radiation the person receives (W/m^2)
% SWRdir_in_top: Direct shortwave radiation the person receives from the top (W/m^2)
% SWRdir_in_bottom: Direct shortwave radiation the person receives from the bottom (W/m^2)
% SWRdir_in_east: Direct shortwave radiation the person receives from the east (W/m^2)
% SWRdir_in_south: Direct shortwave radiation the person receives from the south (W/m^2)
% SWRdir_in_west: Direct shortwave radiation the person receives from the west (W/m^2)
% SWRdir_in_north: Direct shortwave radiation the person receives from the north (W/m^2)
% SWRdiff_Person: Diffuse shortwave radiation the person receives (W/m^2)
% LWR_Person: Longwave radiation the person receives (W/m^2)

MeanRadiantTemperature					=	[];
MeanRadiantTemperature.Tmrt				=	zeros(n,1,m);
MeanRadiantTemperature.BoleanInSun		=	zeros(n,1,m);
MeanRadiantTemperature.SWRdir_Person	=	zeros(n,1,m);
MeanRadiantTemperature.SWRdir_in_top	=	zeros(n,1,m);
MeanRadiantTemperature.SWRdir_in_bottom	=	zeros(n,1,m);
MeanRadiantTemperature.SWRdir_in_east	=	zeros(n,1,m);
MeanRadiantTemperature.SWRdir_in_south	=	zeros(n,1,m);
MeanRadiantTemperature.SWRdir_in_west	=	zeros(n,1,m);
MeanRadiantTemperature.SWRdir_in_north	=	zeros(n,1,m);
MeanRadiantTemperature.SWRdiff_Person	=	zeros(n,1,m);
MeanRadiantTemperature.LWR_Person		=	zeros(n,1,m);

MeanRadiantTemperatureNames	=	{'Tmrt';'BoleanInSun';'SWRdir_Person';'SWRdir_in_top';'SWRdir_in_bottom';...
	'SWRdir_in_east';'SWRdir_in_south';'SWRdir_in_west';'SWRdir_in_north';'SWRdiff_Person';'LWR_Person';};

%--------------------------------------------------------------------------
% Albedo: AlbedoOutput
%--------------------------------------------------------------------------
% TotalUrban: Albedo of the total urban area (-)
% TotalCanyon: Albedo of the total canyon area (-)
% Roof: Albedo of the total roof area (-)

AlbedoOutput			=	[];
AlbedoOutput.TotalUrban	=	zeros(n,1,m);
AlbedoOutput.TotalCanyon=	zeros(n,1,m);
AlbedoOutput.Roof		=	zeros(n,1,m);

%--------------------------------------------------------------------------
%% Outdoor thermal comfort: UTCI (degC)
%--------------------------------------------------------------------------
UTCI	=	zeros(n,1,m);

%--------------------------------------------------------------------------
%% LAI output for varying LAI: LAI_ts
%--------------------------------------------------------------------------
% LAI_R = LAI of roof vegetation (-)
% LAI_G = LAI of ground vegetation (-)
% LAI_T = LAI of tree vegetation (-)

LAI_ts          =	[];
LAI_ts.LAI_R	=	zeros(n,1,m);
LAI_ts.LAI_G	=	zeros(n,1,m);
LAI_ts.LAI_T	=	zeros(n,1,m);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
Opt_Solv = optimoptions('lsqnonlin','Display','off');

for ittm = 1:m

%% Initialize variables
[Gemeotry_m,ParTree,geometry,FractionsRoof,FractionsGround,...
	WallLayers,ParSoilRoof,ParSoilGround,ParInterceptionTree,...
	PropOpticalRoof,PropOpticalGround,PropOpticalWall,PropOpticalTree,...
	ParThermalRoof,ParThermalGround,ParThermalWall,ParThermalTree,...
	ParVegRoof,ParVegGround,ParVegTree,Person]=feval(strcat('data_functions.Data_UEHM_site_',Name_Site),MeteoData,ittm,LAI_TimeSeries);


%% Calculate view factors
[ViewFactor,ViewFactorPoint]=ray_tracing.VFUrbanCanyon(OPTION_RAY,Name_Site,Gemeotry_m,geometry,Person,ParTree);

if FractionsRoof.fimp==1
Vwater.VRoofSoilVeg(1,:,ittm)	=	repmat(0.*ParSoil.Roof.dz,1,1,1);		% Starting with dry soil
Owater.OwRoofSoilVeg(1,:,ittm)	=	0;				% Starting with dry soil
end

if FractionsGround.fimp==1
Vwater.VGroundSoilImp(1,:,ittm)		=	repmat(0.*ParSoil.Ground.dz,1,1,1);	% Starting with dry soil
Owater.OwGroundSoilImp(1,:,ittm)	=	0;				% Starting with dry soil
Vwater.VGroundSoilBare(1,:,ittm)	=	repmat(0.*ParSoil.Ground.dz,1,1,1);	% Starting with dry soil
Owater.OwGroundSoilBare(1,:,ittm)	=	0;				% Starting with dry soil
Vwater.VGroundSoilVeg(1,:,ittm)		=	repmat(0.*ParSoil.Ground.dz,1,1,1);	% Starting with dry soil
Owater.OwGroundSoilVeg(1,:,ittm)	=	0;				% Starting with dry soil
end

% to save initial soil moisture for back-computation
OwaterInitial.OwRoofSoilVeg(1,:,ittm)     = Owater.OwRoofSoilVeg(1,:,ittm);
OwaterInitial.OwGroundSoilImp(1,:,ittm)   = Owater.OwGroundSoilImp(1,:,ittm);
OwaterInitial.OwGroundSoilBare(1,:,ittm)  = Owater.OwGroundSoilBare(1,:,ittm);
OwaterInitial.OwGroundSoilVeg(1,:,ittm)   = Owater.OwGroundSoilVeg(1,:,ittm);
OwaterInitial.OwGroundSoilTot(1,:,ittm)   = Owater.OwGroundSoilTot(1,:,ittm);

for ittn	= 1:n  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SunPosition,MeteoData,HumidityAtm,Anthropogenic,location,ParCalculation]...
	=feval(strcat('data_functions.UEHMForcingData_',Name_SiteFD),MeteoDataRaw,ittn,SoilPotW);

%% Initialize variables
[Gemeotry_m,ParTree,geometry,FractionsRoof,FractionsGround,...
	WallLayers,ParSoilRoof,ParSoilGround,ParInterceptionTree,...
	PropOpticalRoof,PropOpticalGround,PropOpticalWall,PropOpticalTree,...
	ParThermalRoof,ParThermalGround,ParThermalWall,ParThermalTree,...
	ParVegRoof,ParVegGround,ParVegTree,Person]=feval(strcat('data_functions.Data_UEHM_site_',Name_Site),MeteoData,ittm,LAI_TimeSeries);  


	for i=1:size(TempVecNames,1)
		if ittn==1
			TempVec_ittm.(cell2mat(TempVecNames(i)))	=	TempVec.(cell2mat(TempVecNames(i)))(1,:,ittm); 
		else
			TempVec_ittm.(cell2mat(TempVecNames(i)))	=	TempVec.(cell2mat(TempVecNames(i)))(ittn-1,:,ittm);
		end
	end
	
	for i=1:size(IntNames,1)
		if ittn==1
			Int_ittm.(cell2mat(IntNames(i)))	=	Int.(cell2mat(IntNames(i)))(1,:,ittm); 
		else
			Int_ittm.(cell2mat(IntNames(i)))	=	Int.(cell2mat(IntNames(i)))(ittn-1,:,ittm); 
		end
	end
	
	for i=1:size(ExWaterNames,1)
		if ittn==1
			ExWater_ittm.(cell2mat(ExWaterNames(i)))	=	ExWater.(cell2mat(ExWaterNames(i)))(1,:,ittm);
		else
			ExWater_ittm.(cell2mat(ExWaterNames(i)))	=	ExWater.(cell2mat(ExWaterNames(i)))(ittn-1,:,ittm);
		end
	end
	
	VwaterNames	=	fieldnames(Vwater);
	for i=1:size(VwaterNames,1)
		if ittn==1
			Vwater_ittm.(cell2mat(VwaterNames(i)))	=	Vwater.(cell2mat(VwaterNames(i)))(1,:,ittm);
		else
			Vwater_ittm.(cell2mat(VwaterNames(i)))	=	Vwater.(cell2mat(VwaterNames(i)))(ittn-1,:,ittm);
		end
	end
	
	OwaterNames	=	fieldnames(Owater);
	for i=1:size(OwaterNames,1)
		if ittn==1
			Owater_ittm.(cell2mat(OwaterNames(i)))	=	Owater.(cell2mat(OwaterNames(i)))(1,:,ittm);
		else
			Owater_ittm.(cell2mat(OwaterNames(i)))	=	Owater.(cell2mat(OwaterNames(i)))(ittn-1,:,ittm);
		end
	end
	
	for i=1:size(SoilPotWNames,1)
		if ittn==1
			SoilPotW_ittm.(cell2mat(SoilPotWNames(i)))	=	SoilPotW.(cell2mat(SoilPotWNames(i)))(1,:,ittm);
		else
			SoilPotW_ittm.(cell2mat(SoilPotWNames(i)))	=	SoilPotW.(cell2mat(SoilPotWNames(i)))(ittn-1,:,ittm);
		end
	end
	
	for i=1:size(CiCO2LeafNames,1)
		if ittn==1
			CiCO2Leaf_ittm.(cell2mat(CiCO2LeafNames(i)))	=	CiCO2Leaf.(cell2mat(CiCO2LeafNames(i)))(1,:,ittm); 
		else
			CiCO2Leaf_ittm.(cell2mat(CiCO2LeafNames(i)))	=	CiCO2Leaf.(cell2mat(CiCO2LeafNames(i)))(ittn-1,:,ittm); 
		end
	end
	
	for i=1:size(HumidityNames,1)
		if ittn==1
			Humidity_ittm.(cell2mat(HumidityNames(i)))	=	Humidity.(cell2mat(HumidityNames(i)))(1,:,ittm); 
		else
			Humidity_ittm.(cell2mat(HumidityNames(i)))	=	Humidity.(cell2mat(HumidityNames(i)))(ittn-1,:,ittm); 
		end
	end
	
	for i=1:size(TempDampNames,1)
		if ittn==1
			TempDamp_ittm.(cell2mat(TempDampNames(i)))	=	TempDamp.(cell2mat(TempDampNames(i)))(1,:,ittm);
		else
			TempDamp_ittm.(cell2mat(TempDampNames(i)))	=	TempDamp.(cell2mat(TempDampNames(i)))(ittn-1,:,ittm);
		end
	end
	
	for i=1:size(RunonNames,1)
		if ittn==1
			Runon_ittm.(cell2mat(RunonNames(i)))	=	Runon.(cell2mat(RunonNames(i)))(1,:,ittm); 
		else
			Runon_ittm.(cell2mat(RunonNames(i)))	=	Runon.(cell2mat(RunonNames(i)))(ittn-1,:,ittm); 
		end
	end
	
	QinlatNames	=	fieldnames(Qinlat);
	for i=1:size(QinlatNames,1)
		if ittn==1
			Qinlat_ittm.(cell2mat(QinlatNames(i)))	=	Qinlat.(cell2mat(QinlatNames(i)))(1,:,ittm); 
		else
			Qinlat_ittm.(cell2mat(QinlatNames(i)))	=	Qinlat.(cell2mat(QinlatNames(i)))(ittn-1,:,ittm); 
		end
    end
    	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
% Calculate Energy Budget of Roof
	ittn_only	=	1;
	[TR,fvalR,exitflagR]=fSolver_roof(TempVec_ittm,MeteoData,ittn_only,...
	Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,CiCO2Leaf_ittm,Opt_Solv,...
	Gemeotry_m,FractionsRoof,ParSoilRoof,PropOpticalRoof,ParThermalRoof,ParVegRoof,...
	HumidityAtm,Anthropogenic,ParCalculation);

	TempVec.TRoofImp(ittn,1,ittm)		=	TR(1,1);
	TempVec.TRoofVeg(ittn,1,ittm)		=	TR(1,2);
	TempVec.TRoofIntImp(ittn,1,ittm)	=	TR(1,3);
	TempVec.TRoofIntVeg(ittn,1,ittm)	=	TR(1,4);
	
% Calculate Energy Budget of the Canyon
	[TC,fvalC,exitflagC]=fSolver_canyon(TempVec_ittm,Humidity_ittm,MeteoData,ittn_only,...
	Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,CiCO2Leaf_ittm,TempDamp_ittm,ViewFactor,Opt_Solv,...
	Gemeotry_m,ParTree,geometry,FractionsGround,...
	WallLayers,ParSoilGround,ParInterceptionTree,...
	PropOpticalGround,PropOpticalWall,PropOpticalTree,...
	ParThermalGround,ParThermalWall,ParVegGround,ParVegTree,...
	SunPosition,HumidityAtm,Anthropogenic,ParCalculation);

	TempVec.TGroundImp(ittn,1,ittm)		=	TC(1,1);
	TempVec.TGroundBare(ittn,1,ittm)	=	TC(1,2);
	TempVec.TGroundVeg(ittn,1,ittm)		=	TC(1,3);
	TempVec.TTree(ittn,1,ittm)			=	TC(1,6);
	TempVec.TWallSun(ittn,1,ittm)		=	TC(1,4);
	TempVec.TWallShade(ittn,1,ittm)		=	TC(1,5);
	TempVec.TWallIntSun(ittn,1,ittm)	=	TC(1,7);
	TempVec.TWallIntShade(ittn,1,ittm)	=	TC(1,8);
	TempVec.TCanyon(ittn,1,ittm)		=	TC(1,9);
	Humidity.CanyonSpecific(ittn,1,ittm)=	TC(1,10);

% Calculate Energy and Water fluxes Roof
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SWRabsRoofImp,SWRabsRoofVeg,SWRabsTotalRoof,...
		SWRoutRoofImp,SWRoutRoofVeg,SWRoutTotalRoof,...
		SWRinRoofImp,SWRinRoofVeg,SWRinTotalRoof,...
		SWREBRoofImp,SWREBRoofVeg,SWREBTotalRoof,...
		LWRabsRoofVeg,LWRabsRoofImp,LWRabsTotalRoof,...
		LWRoutRoofVeg,LWRoutRoofImp,LWRoutTotalRoof,...
		LWRinRoofImp,LWRinRoofVeg,LWRinTotalRoof,...
		LWREBRoofImp,LWREBRoofVeg,LWREBTotalRoof,...
		HfluxRoofImp,HfluxRoofVeg,HfluxRoof,...
		LEfluxRoofImp,LEfluxRoofVegInt,LEfluxRoofVegPond,...
		LEfluxRoofVegSoil,LTEfluxRoofVeg,LEfluxRoofVeg,LEfluxRoof,...
		G1RoofImp,G2RoofImp,dsRoofImp,G1RoofVeg,G2RoofVeg,dsRoofVeg,G1Roof,G2Roof,dsRoof,...
		raRooftoAtm,rb_LRoof,rap_LRoof,r_soilRoof,rs_sunRoof,rs_shdRoof,...
		EfluxRoofImp,EfluxRoofVegInt,EfluxRoofVegPond,...
		EfluxRoofVegSoil,TEfluxRoofVeg,EfluxRoofVeg,EfluxRoof,...
		... % Water fluxes
		QRoofImp,QRoofVegDrip,QRoofVegPond,LkRoofImp,LkRoofVeg,LkRoof,QRoofVegSoil,RunoffRoofTot,RunonRoofTot,...
		IntRoofImp,IntRoofVegPlant,IntRoofVegGround,dInt_dtRoofImp,dInt_dtRoofVegPlant,dInt_dtRoofVegGround,...
		IntRooftot,dInt_dtRooftot,dVRoofSoilVeg_dt,...
		fRoofVeg,VRoofSoilVeg,OwRoofSoilVeg,OSwRoofSoilVeg,ExWaterRoofVeg_H,SoilPotWRoofVeg_H,SoilPotWRoofVeg_L,ExWaterRoofVeg_L,...
		CiCO2LeafRoofVegSun,CiCO2LeafRoofVegShd,...
		WBRoofVegInVeg,WBRoofVegInGround,WBRoofVegSoil,...
		...
		EBRoofImp,EBRoofVeg,Yroof,WBRoofImp,WBRoofVeg,WBRoofTot]...
		=EB_WB_roof(TR,TempVec_ittm,MeteoData,ittn_only,...
		Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,CiCO2Leaf_ittm,Runon_ittm,...
		Gemeotry_m,FractionsRoof,ParSoilRoof,PropOpticalRoof,ParThermalRoof,ParVegRoof,...
		HumidityAtm,Anthropogenic,ParCalculation);
	
% Calculate Energy and Water fluxes Canyon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SWRin_t,SWRout_t,SWRabs_t,SWRabsDir_t,SWRabsDiff_t,SWREB_t,albedo_canyon,...
		LWRin_t,LWRout_t,LWRabs_t,LWREB_t,...
		HfluxGroundImp,HfluxGroundBare,HfluxGroundVeg,HfluxTree,HfluxGround,...
		EfluxGroundImp,EfluxGroundBarePond,EfluxGroundBareSoil,EfluxGroundVegInt,...
		EfluxGroundVegPond,EfluxGroundVegSoil,TEfluxGroundVeg,EfluxTreeInt,TEfluxTree,...
		EfluxGroundBare,EfluxGroundVeg,EfluxGround,EfluxTree,...
		LEfluxGroundImp,LEfluxGroundBarePond,LEfluxGroundBareSoil,LEfluxGroundVegInt,...
		LEfluxGroundVegPond,LEfluxGroundVegSoil,LTEfluxGroundVeg,LEfluxTreeInt,LTEfluxTree,...
		LEfluxGroundBare,LEfluxGroundVeg,LEfluxGround,LEfluxTree,...
		CiCO2LeafTreeSun,CiCO2LeafTreeShd,CiCO2LeafGroundVegSun,CiCO2LeafGroundVegShd,...
		raCanyontoAtm,rap_can,rap_Htree_In,rb_HGround,rb_LGround,...
		r_soilGroundbare,r_soilGroundveg,alp_soilGroundbare,alp_soilGroundveg,...
		rs_sunGround,rs_shdGround,rs_sunTree,rs_shdTree,...
		Fsun_L,Fshd_L,dw_L,RES_w1,RES_w2,rap_W1_In,rap_W2_In,rap_Zp1,...
		HfluxWallSun,HfluxWallShade,EfluxWallSun,EfluxWallShade,LEfluxWallSun,LEfluxWallShade,HfluxCanyon,LEfluxCanyon,EfluxCanyon,...
		G1WallSun,G2WallSun,dsWallSun,G1WallShade,G2WallShade,dsWallShade,...
		G1GroundImp,TDampGroundImp,G1GroundBare,TDampGroundBare,G1GroundVeg,TDampGroundVeg,GTree,TDampTree,G1Ground,G1Canyon,G2Canyon,...
		dsGroundImp,dsGroundBare,dsGroundVeg,dsTree,dsCanyonAir,Ycanyon,...
		... %Water fluxes
		QTree,IntTree,dInt_dtTree,QGroundVegDrip,IntGroundVegPlant,dInt_dtGroundVegPlant,...
		QGroundImp,IntGroundImp,dInt_dtGroundImp,fGroundImp,QGroundBarePond,IntGroundBare,dInt_dtGroundBare,fGroundBare,...
		QGroundVegPond,IntGroundVegGround,dInt_dtGroundVegGround,fGroundVeg,...
		...
		VGroundSoilImp,OwGroundSoilImp,OSwGroundSoilImp,LkGroundImp,SoilPotWGroundImp_H,SoilPotWGroundImp_L,...
		ExWaterGroundImp_H,ExWaterGroundImp_L,Rd_gimp,TEgveg_imp,TEtree_imp,...
		Egimp_soil,dVGroundSoilImp_dt,Psi_Soil_gimp,Kf_gimp,...
		...
		VGroundSoilBare,OwGroundSoilBare,OSwGroundSoilBare,LkGroundBare,SoilPotWGroundBare_H,SoilPotWGroundBare_L,...
		ExWaterGroundBare_H,ExWaterGroundBare_L,QGroundBareSoil,TEgveg_bare,TEtree_bare,...
		Egbare_Soil,dVGroundSoilBare_dt,Psi_soil_gbare,Kf_gbare,...
		...
		VGroundSoilVeg,OwGroundSoilVeg,OSwGroundSoilVeg,LkGroundVeg,SoilPotWGroundVeg_H,SoilPotWGroundVeg_L,...
		ExWaterGroundVeg_H,ExWaterGroundVeg_L,QGroundVegSoil,TEgveg_veg,TEtree_veg,...
		Egveg_Soil,dVGroundSoilVeg_dt,Psi_soil_gveg,Kf_gveg,...
		...
		Qin_imp,Qin_bare,Qin_veg,Qin_bare2imp,Qin_bare2veg,Qin_imp2bare,Qin_imp2veg,Qin_veg2imp,Qin_veg2bare,...
		...
		VGroundSoilTot,OwGroundSoilTot,OSwGroundSoilTot,LkGround,Rd,dVGroundSoilTot_dt,SoilPotWGroundTot_L,ExWaterGroundTot_L,TEgveg_tot,SoilPotWGroundTot_H,ExWaterGroundTot_H,...
		TEtree_tot,EB_TEtree,EB_TEgveg,WBIndv,WBTot,...
		RunoffGroundTot,RunonGroundTot,Etot,DeepGLk,StorageTot,...
		...
		EBGroundImp,EBGroundBare,EBGroundVeg,EBTree,EBWallSun,EBWallShade,EBWallSunInt,EBWallShadeInt,EBCanyonT,EBCanyonQ,...
		HumidityCan,HumidityAtm,u_Hcan,u_Zref_und,T2m,q2m,e_T2m,RH_T2m,qcan,e_Tcan,RH_Tcan,...
		....
		DHi,Himp_2m,Hbare_2m,Hveg_2m,Hwsun_2m,Hwshade_2m,Hcan_2m,...
		DEi,Eimp_2m,Ebare_soil_2m,Eveg_int_2m,Eveg_soil_2m,TEveg_2m,Ecan_2m]...
		=EB_WB_canyon(TC,TempVec_ittm,MeteoData,ittn_only,...
		Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,...
		CiCO2Leaf_ittm,TempDamp_ittm,Runon_ittm,Qinlat_ittm,ViewFactor,...
		Gemeotry_m,ParTree,geometry,FractionsGround,...
		WallLayers,ParSoilGround,ParInterceptionTree,...
		PropOpticalGround,PropOpticalWall,PropOpticalTree,...
		ParThermalGround,ParThermalWall,ParVegGround,ParVegTree,...
		SunPosition,HumidityAtm,Anthropogenic,ParCalculation);
	
% Outdoor Thermal comfort calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Tmrt,BoleanInSun,SWRdir_Person,SWRdir_in_top,SWRdir_in_bottom,...
	SWRdir_in_east,SWRdir_in_south,SWRdir_in_west,SWRdir_in_north,...
	SWRdiff_Person,LWR_Person]=MRT.MeanRadiantTemperature(SWRout_t,LWRout_t,MeteoData,ViewFactorPoint,...
	ParTree,ParVegTree,geometry,Gemeotry_m,SunPosition,Person);
	
[u_ZPerson]=resistance_functions.WindProfile_PointOutput(Person.HeightWind,...
	Gemeotry_m,ParVegTree,ParTree,MeteoData,FractionsGround,ParVegGround);	

[UTCI_approx]=OTC.UTCI_approx(T2m-273.15,RH_T2m.*100,Tmrt,u_ZPerson);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SWRabs_t.SWRabsTotalUrban	=	geometry.wroof_norm*SWRabsTotalRoof + geometry.wcanyon_norm*SWRabs_t.SWRabsTotalCanyon;
SWRin_t.SWRinTotalUrban		=	geometry.wroof_norm*SWRinTotalRoof + geometry.wcanyon_norm*SWRin_t.SWRinTotalCanyon;
SWRout_t.SWRoutTotalUrban	=	geometry.wroof_norm*SWRoutTotalRoof + geometry.wcanyon_norm*SWRout_t.SWRoutTotalCanyon;
SWREB_t.SWREBTotalUrban		=	geometry.wroof_norm*SWREBTotalRoof + geometry.wcanyon_norm*SWREB_t.SWREBTotalCanyon;

LWRabs_t.LWRabsTotalUrban	=	geometry.wroof_norm*LWRabsTotalRoof + geometry.wcanyon_norm*LWRabs_t.LWRabsTotalCanyon;
LWRin_t.LWRinTotalUrban		=	geometry.wroof_norm*LWRinTotalRoof + geometry.wcanyon_norm*LWRin_t.LWRinTotalCanyon;
LWRout_t.LWRoutTotalUrban	=	geometry.wroof_norm*LWRoutTotalRoof + geometry.wcanyon_norm*LWRout_t.LWRoutTotalCanyon;
LWREB_t.LWREBTotalUrban		=	geometry.wroof_norm*LWREBTotalRoof + geometry.wcanyon_norm*LWREB_t.LWREBTotalCanyon;

HfluxUrban	=	geometry.wroof_norm*HfluxRoof + geometry.wcanyon_norm*HfluxCanyon;
LEfluxUrban	=	geometry.wroof_norm*LEfluxRoof + geometry.wcanyon_norm*LEfluxCanyon;
G1Urban		=	geometry.wroof_norm*G1Roof + geometry.wcanyon_norm*G1Canyon;
G2Urban		=	geometry.wroof_norm*G2Roof + geometry.wcanyon_norm*G2Canyon;
EfluxUrban	=	geometry.wroof_norm*EfluxRoof + geometry.wcanyon_norm*EfluxCanyon;
RunonUrban	=	geometry.wroof_norm*RunonRoofTot + geometry.wcanyon_norm*RunonGroundTot;
RunoffUrban	=	geometry.wroof_norm*RunoffRoofTot + geometry.wcanyon_norm*RunoffGroundTot;
LkUrban		=	geometry.wroof_norm*LkRoof + geometry.wcanyon_norm*LkGround;
	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tracking the success of the energy balance solver
Solver.SuccessRoof(ittn,:,ittm)		=	exitflagR;
Solver.SuccessCanyon(ittn,:,ittm)	=	exitflagC;
Solver.ValuesRoof(ittn,:,ittm)		=	fvalR;
Solver.ValuesCanyon(ittn,:,ittm)	=	fvalC;
Solver.TsolverRoof(ittn,:,ittm)		=	TR;
Solver.TsolverCanyon(ittn,:,ittm)	=	TC;

% Results 2m
Results2m.T2m(ittn,:,ittm)		=	T2m;
Results2m.q2m(ittn,:,ittm)		=	q2m;
Results2m.e_T2m(ittn,:,ittm)	=	e_T2m;
Results2m.RH_T2m(ittn,:,ittm)	=	RH_T2m;
Results2m.qcan(ittn,:,ittm)		=	qcan;
Results2m.e_Tcan(ittn,:,ittm)	=	e_Tcan;
Results2m.RH_Tcan(ittn,:,ittm)	=	RH_Tcan;

% Results energy fluxes 2m
Results2mEnergyFlux.DHi(ittn,:,ittm)			=	DHi;
Results2mEnergyFlux.Himp_2m(ittn,:,ittm)		=	Himp_2m;
Results2mEnergyFlux.Hbare_2m(ittn,:,ittm)		=	Hbare_2m;
Results2mEnergyFlux.Hveg_2m(ittn,:,ittm)		=	Hveg_2m;
Results2mEnergyFlux.Hwsun_2m(ittn,:,ittm)		=	Hwsun_2m;
Results2mEnergyFlux.Hwshade_2m(ittn,:,ittm)		=	Hwshade_2m;
Results2mEnergyFlux.Hcan_2m(ittn,:,ittm)		=	Hcan_2m;

Results2mEnergyFlux.DEi(ittn,:,ittm)			=	DEi;
Results2mEnergyFlux.Eimp_2m(ittn,:,ittm)		=	Eimp_2m;
Results2mEnergyFlux.Ebare_soil_2m(ittn,:,ittm)	=	Ebare_soil_2m;
Results2mEnergyFlux.Eveg_int_2m(ittn,:,ittm)	=	Eveg_int_2m;
Results2mEnergyFlux.Eveg_soil_2m(ittn,:,ittm)	=	Eveg_soil_2m;
Results2mEnergyFlux.TEveg_2m(ittn,:,ittm)		=	TEveg_2m;
Results2mEnergyFlux.Ecan_2m	(ittn,:,ittm)		=	Ecan_2m;


% Anthropogenic time series
Anthropo.Tb(ittn,:,ittm)				=	Anthropogenic.Tb;
Anthropo.Qf_canyon(ittn,:,ittm)			=	Anthropogenic.Qf_canyon;
Anthropo.Qf_roof(ittn,:,ittm)			=	Anthropogenic.Qf_roof;
Anthropo.Waterf_canyonVeg(ittn,:,ittm)	=	Anthropogenic.Waterf_canyonVeg;
Anthropo.Waterf_canyonBare(ittn,:,ittm)	=	Anthropogenic.Waterf_canyonBare;
Anthropo.Waterf_roof(ittn,:,ittm)		=	Anthropogenic.Waterf_roof;

% Albedo
albedo_urban	=	geometry.wcanyon_norm*albedo_canyon + geometry.wroof_norm*PropOpticalRoof.albedo;

AlbedoOutput.TotalUrban(ittn,:,ittm)	=	albedo_urban;
AlbedoOutput.TotalCanyon(ittn,:,ittm)	=	albedo_canyon;
AlbedoOutput.Roof(ittn,:,ittm)			=	PropOpticalRoof.albedo;


% Mean radiant temperature
for i=1:length(MeanRadiantTemperatureNames)
	MeanRadiantTemperature.(cell2mat(MeanRadiantTemperatureNames(i)))(ittn,1,ittm)		=	eval(cell2mat(MeanRadiantTemperatureNames(i)));
end

% Universal thermal comfort index (UTCI)
UTCI(ittn,1,ittm)	=	UTCI_approx;

% Assing LAI for cases with varying LAI
LAI_ts.LAI_R(ittn,1,ittm)	=	ParVegRoof.LAI;
LAI_ts.LAI_G(ittn,1,ittm)	=	ParVegGround.LAI;
LAI_ts.LAI_T(ittn,1,ittm)	=	ParVegTree.LAI;


% Humidity
for i=1:6
	Humidity.(cell2mat(HumidityNames(i)))(ittn,1,ittm)=	HumidityCan.(cell2mat(HumidityNames(i)));
end

for i=7:12
	Humidity.(cell2mat(HumidityNames(i)))(ittn,1,ittm)=	HumidityAtm.(cell2mat(HumidityNames(i)));
end

% Shortwave radiation
for i=1:3
	SWRabs.(cell2mat(SWRabsNames(i)))(ittn,1,ittm)	=	eval(cell2mat(SWRabsNames(i)));
	SWRin.(cell2mat(SWRinNames(i)))(ittn,1,ittm)	=	eval(cell2mat(SWRinNames(i)));
	SWRout.(cell2mat(SWRoutNames(i)))(ittn,1,ittm)	=	eval(cell2mat(SWRoutNames(i)));
	SWREB.(cell2mat(SWREBNames(i)))(ittn,1,ittm)	=	eval(cell2mat(SWREBNames(i)));
end
for i=4:size(SWRabsNames,1)
	SWRabs.(cell2mat(SWRabsNames(i)))(ittn,1,ittm)	=	SWRabs_t.(cell2mat(SWRabsNames(i)))(1,1);
	SWRin.(cell2mat(SWRinNames(i)))(ittn,1,ittm)	=	SWRin_t.(cell2mat(SWRinNames(i)))(1,1);
	SWRout.(cell2mat(SWRoutNames(i)))(ittn,1,ittm)	=	SWRout_t.(cell2mat(SWRoutNames(i)))(1,1);
	SWREB.(cell2mat(SWREBNames(i)))(ittn,1,ittm)	=	SWREB_t.(cell2mat(SWREBNames(i)))(1,1);
end

% Longwave radiation
for i=1:3
	LWRabs.(cell2mat(LWRabsNames(i)))(ittn,1,ittm)	=	eval(cell2mat(LWRabsNames(i)));
	LWRin.(cell2mat(LWRinNames(i)))(ittn,1,ittm)	=	eval(cell2mat(LWRinNames(i)));
	LWRout.(cell2mat(LWRoutNames(i)))(ittn,1,ittm)	=	eval(cell2mat(LWRoutNames(i)));
	LWREB.(cell2mat(LWREBNames(i)))(ittn,1,ittm)	=	eval(cell2mat(LWREBNames(i)));
end
for i=4:size(LWRabsNames,1)
	LWRabs.(cell2mat(LWRabsNames(i)))(ittn,1,ittm)	=	LWRabs_t.(cell2mat(LWRabsNames(i)))(1,1);
	LWRin.(cell2mat(LWRinNames(i)))(ittn,1,ittm)	=	LWRin_t.(cell2mat(LWRinNames(i)))(1,1);
	LWRout.(cell2mat(LWRoutNames(i)))(ittn,1,ittm)	=	LWRout_t.(cell2mat(LWRoutNames(i)))(1,1);
	LWREB.(cell2mat(LWREBNames(i)))(ittn,1,ittm)	=	LWREB_t.(cell2mat(LWREBNames(i)))(1,1);
end

% Sensible heat
for i=1:size(HfluxNames,1)
	Hflux.(cell2mat(HfluxNames(i)))(ittn,1,ittm)	=	eval(cell2mat(HfluxNames(i)));
end

% Latent heat
for i=1:size(LEfluxNames,1)
	LEflux.(cell2mat(LEfluxNames(i)))(ittn,1,ittm)	=	eval(cell2mat(LEfluxNames(i)));
end

% Ground heat flux
for i=1:size(GfluxNames,1)
	Gflux.(cell2mat(GfluxNames(i)))(ittn,1,ittm)	=	eval(cell2mat(GfluxNames(i)));
end

% Heat Storage
for i=1:size(dStorageNames,1)
	dStorage.(cell2mat(dStorageNames(i)))(ittn,1,ittm)	=	eval(cell2mat(dStorageNames(i)));
end

% Dampening temperature
for i=1:size(TempDampNames,1)
	TempDamp.(cell2mat(TempDampNames(i)))(ittn,1,ittm)	=	eval(cell2mat(TempDampNames(i)));
end

for i=1:size(RESNames,1)
	RES.(cell2mat(RESNames(i)))(ittn,1,ittm)	=	eval(cell2mat(RESNames(i)));
end

for i=1:size(EfluxNames,1)
	Eflux.(cell2mat(EfluxNames(i)))(ittn,1,ittm)	=	eval(cell2mat(EfluxNames(i)));
end

for i=1:size(RunoffNames,1)
	Runoff.(cell2mat(RunoffNames(i)))(ittn,1,ittm)	=	eval(cell2mat(RunoffNames(i)));
end

for i=1:size(RunonNames,1)
	Runon.(cell2mat(RunonNames(i)))(ittn,1,ittm)	=	eval(cell2mat(RunonNames(i)));
end

for i=1:size(LeakageNames,1)
	Leakage.(cell2mat(LeakageNames(i)))(ittn,1,ittm)	=	eval(cell2mat(LeakageNames(i)));
end

for i=1:size(IntNames,1)
	Int.(cell2mat(IntNames(i)))(ittn,1,ittm)	=	eval(cell2mat(IntNames(i)));
end

for i=1:size(dInt_dtNames,1)
	dInt_dt.(cell2mat(dInt_dtNames(i)))(ittn,1,ittm)	=	eval(cell2mat(dInt_dtNames(i)));
end

for i=1:size(InfiltrationNames,1)
	Infiltration.(cell2mat(InfiltrationNames(i)))(ittn,1,ittm)=	eval(cell2mat(InfiltrationNames(i)));
end

Vwater.VRoofSoilVeg(ittn,:,ittm)	=	VRoofSoilVeg;
Vwater.VGroundSoilImp(ittn,:,ittm)	=	VGroundSoilImp;
Vwater.VGroundSoilBare(ittn,:,ittm)	=	VGroundSoilBare;
Vwater.VGroundSoilVeg(ittn,:,ittm)	=	VGroundSoilVeg;
Vwater.VGroundSoilTot(ittn,:,ittm)	=	VGroundSoilTot;

dVwater_dt.dVRoofSoilVeg_dt(ittn,:,ittm)	=	dVRoofSoilVeg_dt;
dVwater_dt.dVGroundSoilImp_dt(ittn,:,ittm)	=	dVGroundSoilImp_dt;
dVwater_dt.dVGroundSoilBare_dt(ittn,:,ittm)	=	dVGroundSoilBare_dt;
dVwater_dt.dVGroundSoilVeg_dt(ittn,:,ittm)	=	dVGroundSoilVeg_dt;
dVwater_dt.dVGroundSoilTot_dt(ittn,:,ittm)	=	dVGroundSoilTot_dt;

% OwaterAfterExtraction.OwRoofSoilVeg(ittn,:,ittm)      =	OwRoofSoilVeg;
% OwaterAfterExtraction.OwGroundSoilImp(ittn,:,ittm)	=	OwGroundSoilImp;
% OwaterAfterExtraction.OwGroundSoilBare(ittn,:,ittm)   =	OwGroundSoilBare;
% OwaterAfterExtraction.OwGroundSoilVeg(ittn,:,ittm)	=	OwGroundSoilVeg;
% OwaterAfterExtraction.OwGroundSoilTot(ittn,:,ittm)	=	OwGroundSoilTot;

Owater.OwRoofSoilVeg(ittn,:,ittm)	=	OwRoofSoilVeg;
Owater.OwGroundSoilImp(ittn,:,ittm)	=	OwGroundSoilImp;
Owater.OwGroundSoilBare(ittn,:,ittm)=	OwGroundSoilBare;
Owater.OwGroundSoilVeg(ittn,:,ittm)	=	OwGroundSoilVeg;
Owater.OwGroundSoilTot(ittn,:,ittm)	=	OwGroundSoilTot;

% Fixed soil moisture
if ParSoilRoof.FixSM_R==1
    ReplaceVal_R    =   ParSoil.Roof.O33;
    SMReplace_R     =   false(ParSoilRoof.ms,1);
    SMReplace_R(ParSoilRoof.FixSM_LayerStart_R:ParSoilRoof.FixSM_LayerEnd_R,1) = true;
    
    Owater.OwRoofSoilVeg(ittn,(SMReplace_R & OwRoofSoilVeg'<ReplaceVal_R),ittm)	=	ReplaceVal_R;
end

if ParSoilGround.FixSM_G==1
    ReplaceVal_G    =   ParSoil.Ground.O33;
    SMReplace_G     =   false(ParSoilGround.ms,1);
    SMReplace_G(ParSoilGround.FixSM_LayerStart_G:ParSoilGround.FixSM_LayerEnd_G,1) = true;
 
    Owater.OwGroundSoilImp(ittn,(SMReplace_G & OwGroundSoilImp'<ReplaceVal_G),ittm)	=	ReplaceVal_G;	
    Owater.OwGroundSoilBare(ittn,(SMReplace_G & OwGroundSoilBare'<ReplaceVal_G),ittm)=	ReplaceVal_G;	
    Owater.OwGroundSoilVeg(ittn,(SMReplace_G & OwGroundSoilVeg'<ReplaceVal_G),ittm)	=	ReplaceVal_G;	
    Owater.OwGroundSoilTot(ittn,(SMReplace_G & OwGroundSoilTot'<ReplaceVal_G),ittm)	=	ReplaceVal_G;	
end

OSwater.OSwRoofSoilVeg(ittn,:,ittm)	=	OSwRoofSoilVeg;
OSwater.OSwGroundSoilImp(ittn,:,ittm)=	OSwGroundSoilImp;
OSwater.OSwGroundSoilBare(ittn,:,ittm)=	OSwGroundSoilBare;
OSwater.OSwGroundSoilVeg(ittn,:,ittm)=	OSwGroundSoilVeg;
OSwater.OSwGroundSoilTot(ittn,:,ittm)=	OSwGroundSoilTot;

% Lateral soil water flux
Qinlat.Qin_bare2imp(ittn,:,ittm)=	Qin_bare2imp;
Qinlat.Qin_veg2imp(ittn,:,ittm)	=	Qin_veg2imp;
Qinlat.Qin_veg2bare(ittn,:,ittm)=	Qin_veg2bare;
Qinlat.Qin_imp2bare(ittn,:,ittm)=	Qin_imp2bare;
Qinlat.Qin_bare2veg(ittn,:,ittm)=	Qin_bare2veg;
Qinlat.Qin_imp2veg(ittn,:,ittm)	=	Qin_imp2veg;
Qinlat.Qin_imp(ittn,:,ittm)		=	Qin_imp;
Qinlat.Qin_bare(ittn,:,ittm)	=	Qin_bare;
Qinlat.Qin_veg(ittn,:,ittm)		=	Qin_veg;


for i=1:size(ExWaterNames,1)
	ExWater.(cell2mat(ExWaterNames(i)))(ittn,:,ittm)=	eval(cell2mat(ExWaterNames(i)));
end

for i=1:size(SoilPotWNames,1)
	SoilPotW.(cell2mat(SoilPotWNames(i)))(ittn,1,ittm)=	eval(cell2mat(SoilPotWNames(i)));
end

for i=1:size(CiCO2LeafNames,1)
	CiCO2Leaf.(cell2mat(CiCO2LeafNames(i)))(ittn,1,ittm)=	eval(cell2mat(CiCO2LeafNames(i)));
end

for i=1:size(EBNames,1)
	EB.(cell2mat(EBNames(i)))(ittn,1,ittm)=	eval(cell2mat(EBNames(i)));
end

WBRoofVegSoil	=	sum(WBRoofVegSoil);
for i=1:size(WBRoofNames,1)
	WBRoof.(cell2mat(WBRoofNames(i)))(ittn,:,ittm)=	eval(cell2mat(WBRoofNames(i)));
end

WBCanyonIndv.WB_In_tree(ittn,:,ittm)	=	WBIndv.WB_In_tree;
WBCanyonIndv.WB_In_gveg(ittn,:,ittm)	=	WBIndv.WB_In_gveg;
WBCanyonIndv.WB_In_gimp(ittn,:,ittm)	=	WBIndv.WB_In_gimp;
WBCanyonIndv.WB_In_gbare(ittn,:,ittm)	=	WBIndv.WB_In_gbare;
WBCanyonIndv.WB_Pond_gveg(ittn,:,ittm)	=	WBIndv.WB_Pond_gveg;
WBCanyonIndv.WB_Soil_gimp(ittn,:,ittm)	=	nansum(WBIndv.WB_Soil_gimp);
WBCanyonIndv.WB_Soil_gbare(ittn,:,ittm)	=	sum(WBIndv.WB_Soil_gbare);
WBCanyonIndv.WB_Soil_gveg(ittn,:,ittm)	=	sum(WBIndv.WB_Soil_gveg);

WBCanyonTot.WBsurf_tree(ittn,:,ittm)	=	WBTot.WBsurf_tree;
WBCanyonTot.WBsurf_imp(ittn,:,ittm)		=	WBTot.WBsurf_imp;
WBCanyonTot.WBsurf_bare(ittn,:,ittm)	=	WBTot.WBsurf_bare;
WBCanyonTot.WBsurf_veg(ittn,:,ittm)		=	WBTot.WBsurf_veg;
WBCanyonTot.WBsoil_imp(ittn,:,ittm)		=	WBTot.WBsoil_imp;
WBCanyonTot.WBsoil_bare(ittn,:,ittm)	=	WBTot.WBsoil_bare;
WBCanyonTot.WBsoil_veg(ittn,:,ittm)		=	WBTot.WBsoil_veg;
WBCanyonTot.WBimp_tot(ittn,:,ittm)		=	WBTot.WBimp_tot;
WBCanyonTot.WBbare_tot(ittn,:,ittm)		=	WBTot.WBbare_tot;
WBCanyonTot.WBveg_tot(ittn,:,ittm)		=	WBTot.WBveg_tot;
WBCanyonTot.WBcanyon_flux(ittn,:,ittm)	=	WBTot.WBcanyon_flux;
WBCanyonTot.WBtree_level(ittn,:,ittm)	=	WBTot.WBtree_level;
WBCanyonTot.WBground_level(ittn,:,ittm)	=	WBTot.WBground_level;
WBCanyonTot.WBsoil_level(ittn,:,ittm)	=	WBTot.WBsoil_level;
WBCanyonTot.WBcanyon_level(ittn,:,ittm)	=	WBTot.WBcanyon_level;

for i=1:size(EBNames,1)
	EB.(cell2mat(EBNames(i)))(ittn,1,ittm)=	eval(cell2mat(EBNames(i)));
end

for i=1:size(WindNames,1)
	Wind.(cell2mat(WindNames(i)))(ittn,1,ittm)=	eval(cell2mat(WindNames(i)));
end


if mod(ittn,10)==0
disp(strcat('iter=',num2str(ittn),' iterm=',num2str(ittm)));
end 

end

Gemeotry_m_Out(ittm)			=	Gemeotry_m;
ParTree_Out(ittm)				=	ParTree;
geometry_Out(ittm)				=	geometry;
FractionsRoof_Out(ittm)			=	FractionsRoof;
FractionsGround_Out(ittm)		=	FractionsGround;
WallLayers_Out(ittm)			=	WallLayers;
ParSoilRoof_Out(ittm)			=	ParSoilRoof;
ParSoilGround_Out(ittm)			=	ParSoilGround;
ParInterceptionTree_Out(ittm)	=	ParInterceptionTree;
PropOpticalRoof_Out(ittm)		=	PropOpticalRoof;
PropOpticalGround_Out(ittm)		=	PropOpticalGround;
PropOpticalWall_Out(ittm)		=	PropOpticalWall;
PropOpticalTree_Out(ittm)		=	PropOpticalTree;
ParThermalRoof_Out(ittm)		=	ParThermalRoof;
ParThermalGround_Out(ittm)		=	ParThermalGround;
ParThermalWall_Out(ittm)		=	ParThermalWall;
ParThermalTree_Out(ittm)		=	ParThermalTree;
ParVegRoof_Out(ittm)			=	ParVegRoof;
ParVegGround_Out(ittm)			=	ParVegGround;
ParVegTree_Out(ittm)			=	ParVegTree;
ParCalculation_Out(ittm)		=	ParCalculation;


end

Computational_Time =toc;
%profile off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('COMPUTATIONAL TIME [s] ')
disp(Computational_Time)
disp(' COMPUTATIONAL TIME [ms/cycle] ')
disp(1000*Computational_Time/ittn)


% Because we overwrite the initial temperature at first time step we get a
% wrong temperature storage and hence a wrong energy balance in the first
% time step.
for i=1:size(EBNames,1)
	EB.(cell2mat(EBNames(i)))(1,1,:)	=	0;
end

Zatm = MeteoData.Zatm;

% Plot and calculate radiation and energy balance
[WaterFluxRoof,WaterFluxCan,WaterFluxUrban]=WaterBalanceComponents(MeteoDataRaw,...
    Runon,Leakage,LEflux,dVwater_dt,OwaterInitial,Owater,dInt_dt,Int,Anthropo,...
    ParSoil,ParCalculation_Out,FractionsRoof_Out,FractionsGround_Out,geometry_Out,1);

UrbanClimateVariables(TempVec,UTCI,Results2m,MeteoDataRaw,MeanRadiantTemperature,...
    FractionsGround_Out,FractionsRoof_Out,ParTree_Out,1);

[EnergyFluxUrban,EnergyFluxCan,EnergyFluxRoof]=PlanAreaEnergyBalanceCalculation(ViewFactor,MeteoDataRaw,...
    SWRin,SWRout,SWRabs,LWRin,LWRout,LWRabs,LEflux,Hflux,Gflux,...
    geometry_Out,FractionsGround_Out,PropOpticalRoof_Out,Anthropo,1);



save(['Calculation',NameOutput],'Solver','TempVec','Humidity','SWRabs','SWRin','SWRout','SWREB','LWRabs','LWRin','LWRout',...
	'LWREB','Hflux','LEflux','Gflux','dStorage','RES','Eflux','Runoff','Runon','Leakage',...
	'Int','dInt_dt','Infiltration','Vwater','dVwater_dt','Owater',...
    'OwaterInitial','OSwater','ExWater','SoilPotW',...
	'CiCO2Leaf','WBRoof','WBCanyonIndv','WBCanyonTot','EB','Wind','TempDamp','Qinlat','Results2m',...
	'n','m','Name_Site','MeteoDataRaw','Anthropo',...
	'Gemeotry_m_Out','ParTree_Out','geometry_Out','FractionsRoof_Out','FractionsGround_Out',...
	'WallLayers_Out','ParSoilRoof_Out','ParSoilGround_Out','ParInterceptionTree_Out',...
	'PropOpticalRoof_Out','PropOpticalGround_Out','PropOpticalWall_Out','PropOpticalTree_Out',...
	'ParThermalRoof_Out','ParThermalGround_Out','ParThermalWall_Out','ParThermalTree_Out','ParCalculation_Out',...
	'ParVegRoof_Out','ParVegGround_Out','ParVegTree_Out','LAI_ts','Results2mEnergyFlux','MeanRadiantTemperature','Zatm','UTCI',...
    'AlbedoOutput','ViewFactor','EnergyFluxUrban','EnergyFluxCan','EnergyFluxRoof',...
    'WaterFluxRoof','WaterFluxCan','WaterFluxUrban','ParSoil')


%% Check the energy balance and calculation
EnergyBalanceCheck




