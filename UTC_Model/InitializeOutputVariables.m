function[TempVec,TempVecNames,TempDamp,TempDampNames,Humidity,HumidityNames,...
    SWRabs,SWRabsNames,SWRin,SWRinNames,SWRout,SWRoutNames,SWREBNames,SWREB,...
    LWRabs,LWRabsNames,LWRin,LWRinNames,LWRout,LWRoutNames,LWREB,LWREBNames,...
    Hflux,HfluxNames,LEflux,LEfluxNames,Gflux,GfluxNames,dStorage,dStorageNames,...
    RES,RESNames,Eflux,EfluxNames,Runoff,RunoffNames,Runon,RunonNames,Leakage,LeakageNames,...
    Int,IntNames,dInt_dt,dInt_dtNames,Infiltration,InfiltrationNames,...
    Vwater,dVwater_dt,Owater,OSwater,Qinlat,QinlatNames,ExWater,ExWaterNames,SoilPotW,SoilPotWNames,...
    CiCO2Leaf,CiCO2LeafNames,WBRoof,WBRoofNames,WBCanyonIndv,WBCanyonIndvNames,...
    WBCanyonTot,WBCanyonTotNames,EB,EBNames,Wind,WindNames,Solver,Results2m,Results2mNames,...
    Results2mEnergyFlux,Results2mEnergyFluxNames,MeanRadiantTemperature,MeanRadiantTemperatureNames,...
    AlbedoOutput,AlbedoOutputNames,UTCI,LAI_ts,LAI_tsNames,TempVecB,TempVecBNames,HbuildInt,HbuildIntNames,...
    LEbuildInt,LEbuildIntNames,GbuildInt,GbuildIntNames,SWRabsB,SWRabsBNames,...
    LWRabsB,LWRabsBNames,BEMWasteHeat,BEMWasteHeatNames,BEMEnergyUse,BEMEnergyUseNames,...
    HumidityBuilding,HumidityBuildingNames,ParACHeat,ParACHeatNames,...
    ...
    MeteoData,ParSoil]=...
    InitializeOutputVariables(MeteoDataRaw,n,m,Name_Site,Name_SiteFD,LAI_TimeSeries)


% Initialize the vectors in which the simulated variables will be saved in
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializing vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% Temperature: TempVec
%--------------------------------------------------------------------------
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
	TempVec.(cell2mat(TempVecNames(i)))(1,:,:)	=	MeteoData.Tatm;
end

TempVec.Tatm	=	repmat(MeteoDataRaw.T_atm(:,1),1,1,m); % Temperature atmosphere(measured)


% Dampening temperature: TempDamp
%--------------------------------------------------------------------------
% TGroundImp	=	Dampening temperature ground impervious area [K]
% TGroundBare	=	Dampening temperature ground bare area [K]
% TGroundVeg	=	Dampening temperature ground vegetated area [K]
% TTree			=	Dampening temperature tree canopy [K]
% TDampGroundBuild =	Dampening temperature of ground in building interior [K]

TempDampNames	=	{'TDampGroundImp';'TDampGroundBare';'TDampGroundVeg';'TDampTree';'TDampGroundBuild'};

for i=1:size(TempDampNames,1)
	TempDamp.(cell2mat(TempDampNames(i)))		=	zeros(n,1,m);
	TempDamp.(cell2mat(TempDampNames(i)))(1,:,:)=	mean(MeteoDataRaw.T_atm,"omitnan");
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
					'SWRabsTotalGround';'SWRabsTotalCanyon';'SWRabsTotalUrban';...
                    'SWRabsWallSunExt'; 'SWRabsWallShadeExt';...
                    'SWRabsWallSunTransmitted';'SWRabsWallShadeTransmitted'};

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
end

for i=1:size(SWRinNames,1)
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
					'HfluxCanyon';'HfluxUrban';'dS_H_air'};

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
					'LEfluxWallShade';'LEfluxCanyon';'LEfluxUrban';'dS_LE_air'};



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
	
RESNames	=	{'raRooftoAtm';'raCanyontoAtmOrig';'rap_LRoof';'rb_LRoof';'r_soilRoof';'rs_sunRoof';'rs_shdRoof';...
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
QinlatNames	=	{'Qin_bare2imp';'Qin_veg2imp';'Qin_veg2bare';'Qin_imp2bare';...
					'Qin_bare2veg';'Qin_imp2veg';'Qin_imp';'Qin_bare';'Qin_veg'};

for i=1:size(QinlatNames,1)
	Qinlat.(cell2mat(QinlatNames(i)))	=	zeros(n,ParSoil.Ground.ms,m);
end


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
% Success: Bolean indicating convergence of solution of energy balance
% ValuesEB: Energy balance closure for the different equations (W/m^2)
% Tsolver: Temperatures and humidity of different canyon factes and air (K), (kg/kg)

Solver				=	[];
Solver.Success	=	zeros(n,1,m);
Solver.ValuesEB	=	zeros(n,22,m);
Solver.Tsolver	=	zeros(n,22,m);
Solver.YfunctionOutput	=	zeros(n,22,m);

%--------------------------------------------------------------------------
%% Temperature and humidity at 2m canyon height: Results2m
%--------------------------------------------------------------------------
% T2m: 2m air temperature (K)
% q2m: 2m specific humidity (kg/kg)
% e_T2m: 2m vapor pressure (Pa)
% RH_T2m: 2m relative humidity (-)


Results2mNames	=	{'T2m';'q2m';'e_T2m';'RH_T2m';'qcan';'e_Tcan';'RH_Tcan'};

for i=1:size(Results2mNames,1)
	Results2m.(cell2mat(Results2mNames(i)))	=	zeros(n,1,m);
end

Results2m.T2m(1,:,:) =	MeteoData.Tatm;
Results2m.q2m(1,:,:) =	MeteoData.q_atm;

%--------------------------------------------------------------------------
%% Energy fluxes at 2m canyon height: Results2mEnergyFlux
%--------------------------------------------------------------------------
Results2mEnergyFluxNames	=	{'DHi';'Himp_2m';'Hbare_2m';'Hveg_2m';'Hwsun_2m';'Hwshade_2m';'Hcan_2m';...
    'DEi';'Eimp_2m';'Ebare_soil_2m';'Eveg_int_2m';'Eveg_soil_2m';'TEveg_2m';'Ecan_2m'};

for i=1:size(Results2mEnergyFluxNames,1)
	Results2mEnergyFlux.(cell2mat(Results2mEnergyFluxNames(i)))	=	zeros(n,1,m);
end


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

MeanRadiantTemperatureNames	=	{'Tmrt';'BoleanInSun';'SWRdir_Person';'SWRdir_in_top';'SWRdir_in_bottom';...
	'SWRdir_in_east';'SWRdir_in_south';'SWRdir_in_west';'SWRdir_in_north';'SWRdiff_Person';'LWR_Person'};

for i=1:size(MeanRadiantTemperatureNames,1)
	MeanRadiantTemperature.(cell2mat(MeanRadiantTemperatureNames(i)))	=	zeros(n,1,m);
end


%--------------------------------------------------------------------------
% Albedo: AlbedoOutput
%--------------------------------------------------------------------------
% TotalUrban: Albedo of the total urban area (-)
% TotalCanyon: Albedo of the total canyon area (-)
% Roof: Albedo of the total roof area (-)

AlbedoOutputNames	=	{'TotalUrban';'TotalCanyon';'Roof'};

for i=1:size(AlbedoOutputNames,1)
	AlbedoOutput.(cell2mat(AlbedoOutputNames(i)))	=	zeros(n,1,m);
end


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

LAI_tsNames	=	{'LAI_R';'LAI_G';'LAI_T'};

for i=1:size(LAI_tsNames,1)
	LAI_ts.(cell2mat(LAI_tsNames(i)))	=	zeros(n,1,m);
end


%% Initialize variables for building energy model
%--------------------------------------------------------------------------
% Initialize temperature calculation vector for building interiors
TempVecBNames	=	{'Tceiling';'Tinwallsun';'Tinwallshd';'Twindows';'Tinground';'Tintmass';'Tbin';'qbin'};

for i=1:size(TempVecBNames,1)
	TempVecB.(cell2mat(TempVecBNames(i)))			=	zeros(n,1,m);
	TempVecB.(cell2mat(TempVecBNames(i)))(1,:,:)	=	MeteoData.Tatm;
end

TempVecB.qbin(1,:,:) = HumidityAtm.AtmSpecific;


% Initialize energy flux outputs for building interior
% Sensible heat
HbuildIntNames	=	{'HBinRoof';'HbinWallSun';'HbinWallshd';'HBinGround';'HbinIntMass';...
    'HbuildInSurf';'Hvent';'Hequip';'Hpeople';'H_AC_Heat';'dSH_air'};

for i=1:size(HbuildIntNames,1)
	HbuildInt.(cell2mat(HbuildIntNames(i)))			=	zeros(n,1,m);
end

% Latent heat
LEbuildIntNames	=	{'LEvent';'LEequip';'LEpeople';'LE_AC_Heat';'dSLE_air'};

for i=1:size(LEbuildIntNames,1)
	LEbuildInt.(cell2mat(LEbuildIntNames(i)))			=	zeros(n,1,m);
end

% Conductive heat fluxes
GbuildIntNames	=	{'G2Roof';'G2WallSun';'G2WallShade';'Gfloor';'dSinternalMass'};

for i=1:size(GbuildIntNames,1)
	GbuildInt.(cell2mat(GbuildIntNames(i)))			=	zeros(n,1,m);
end

% Absorbed shortave radiation
SWRabsBNames	=	{'SWRabsCeiling';'SWRabsWallsun';'SWRabsWallshd';'SWRabsGround';'SWRabsInternalMass'};

for i=1:size(SWRabsBNames,1)
	SWRabsB.(cell2mat(SWRabsBNames(i)))			=	zeros(n,1,m);
end

% Absorbed longwave radiation
LWRabsBNames	=	{'LWRabsCeiling';'LWRabsWallsun';'LWRabsWallshd';'LWRabsGround';'LWRabsInternalMass'};

for i=1:size(LWRabsBNames,1)
	LWRabsB.(cell2mat(LWRabsBNames(i)))			=	zeros(n,1,m);
end

% Wasterheat from building interiors emitted into canyon
BEMWasteHeatNames	=	{'SensibleFromAC_Can';'LatentFromAC_Can';'WaterFromAC_Can';...
    'SensibleFromHeat_Can';'LatentFromHeat_Can';'SensibleFromVent_Can';'LatentFromVent_Can';'TotAnthInput_URB'};

for i=1:size(BEMWasteHeatNames,1)
	BEMWasteHeat.(cell2mat(BEMWasteHeatNames(i))) = zeros(n,1,m);
end

% Building energy use, Wh -> energy used (Wh) in the calculation time step.
% I.e. 1h run time = 400 Wh for cooling and heating in that hour, but if
% timestep is 0.5h, it will be 200 Wh. Watt hour is the unit of energy used
% and W is the power. Energy = Time * Power
BEMEnergyUseNames	=	{'EnergyForAC';'EnergyForAC_H';'EnergyForAC_LE';'EnergyForHeating'};

for i=1:size(BEMEnergyUseNames,1)
	BEMEnergyUse.(cell2mat(BEMEnergyUseNames(i))) = zeros(n,1,m);
end

% Humidity within building
HumidityBuildingNames	=	{'qbin';'esatbin';'ebin';'RHbin'};

for i=1:size(HumidityBuildingNames,1)
	HumidityBuilding.(cell2mat(HumidityBuildingNames(i))) = zeros(n,1,m);
end

% Time varying AC parameters
ParACHeatNames	=	{'AC_on';'AC_onCool';'AC_onDehum';'Heat_on'};

for i=1:size(ParACHeatNames,1)
	ParACHeat.(cell2mat(ParACHeatNames(i))) = zeros(n,1,m);
end







