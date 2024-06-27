%%%%%%%%%% RUN TIME SERIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile('+data_functions', 'TMYRiyadh_RadPart.mat'))

% Decide if a varying LAI timeseries is provided [1] or not [0]
[LAI_TimeSeries]=data_functions.VaryingLAIInput(0,'LAI_Zurich_Area'); 

n			=	size(TMYRiyadh,1); % Calculation length, there is no need to change this
m			=	1;					% Length for sensitivity analysis
Name_Site	=	'RY_LCZ3';	% Name for Data_UEHM_site
Name_SiteFD	=	'RY_LCZ3';		% Name for UEHMForcingData
Name_ViewFactor = 'RY_LCZ3';
OPTION_RAY	=	1; % Load precalculated view factors [1], Recalculate view factors [0]

NameOutput	=	'RY_LCZ3';

% Specify which outputs should be saved
% 1 = Essential energy flux and climate output
% 2 = Extended energy flux and climate outputs
% 3 = Extended outputs
OutputsToSave = 3; 


%% Meteo data
% LWR_in [W/m2], SAB1_in [W/m2], SAB2_in [W/m2], SAD1_in [W/m2], SAD2_in [W/m2]
% T_atm	[K], windspeed_u[m/s, pressure_atm [Pa], rain [mm/h], rel_humidity [-]
TMYRiyadh.WindSpeedms(TMYRiyadh.WindSpeedms(1:n,:)==0) = 0.01;	% Wind speed cannot be 0 otherwise the resistance function fails

MeteoDataRaw	=	struct('LWR_in',TMYRiyadh.LWRin(1:n,:),'SAB1_in',TMYRiyadh.SAB1(1:n,:),...
					'SAB2_in',TMYRiyadh.SAB2(1:n,:),'SAD1_in',TMYRiyadh.SAD1(1:n,:),...
					'SAD2_in',TMYRiyadh.SAD2(1:n,:),'T_atm',TMYRiyadh.TairC(1:n,:)+273.15,...
					'windspeed_u',TMYRiyadh.WindSpeedms(1:n,:),'pressure_atm',TMYRiyadh.PressurePa(1:n,:),...
					'rain',TMYRiyadh.Rain(1:n,:),'rel_humidity',...
					TMYRiyadh.RH(1:n,:)./100,'Date',TMYRiyadh.Time(1:n,:));


%% Calculation starts here. No need to change anything after this point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
clearvars -except MeteoDataRaw n m Name_Site OPTION_RAY NameOutput Name_SiteFD LAI_TimeSeries Name_ViewFactor OutputsToSave

RESPreCalc      = 1; % Pre-calculate stomatal resistance for faster computation (based on temperature of previous time step)
fconvPreCalc    = 0; % Pre-calculate convection efficiency enhancement factor for faster computation (based on temperature of previous time step)

% Intialize output variables
[TempVec,TempVecNames,TempDamp,TempDampNames,Humidity,HumidityNames,...
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
    HumidityBuilding,HumidityBuildingNames,ParACHeat_ts,ParACHeatNames,...
    ...
    MeteoData,ParSoil]=...
    InitializeOutputVariables(MeteoDataRaw,n,m,Name_Site,Name_SiteFD,LAI_TimeSeries);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic


for ittm = 1:m

%% Initialize variables
[Gemeotry_m,ParTree,geometry,FractionsRoof,FractionsGround,...
	WallLayers,ParSoilRoof,ParSoilGround,ParInterceptionTree,...
	PropOpticalRoof,PropOpticalGround,PropOpticalWall,PropOpticalTree,...
	ParThermalRoof,ParThermalGround,ParThermalWall,ParThermalTree,...
	ParVegRoof,ParVegGround,ParVegTree,Person,...
    PropOpticalIndoors,ParThermalBulidingInt,ParWindows,ParHVAC,BEM_on]...
    =feval(strcat('data_functions.Data_UEHM_site_',Name_Site),MeteoData,ittm,LAI_TimeSeries);


%% Calculate view factors
[ViewFactor,ViewFactorPoint]=ray_tracing.VFUrbanCanyon(OPTION_RAY,Name_ViewFactor,Gemeotry_m,geometry,Person,ParTree);

% Starting with dry soil if roof is fully impervious
if FractionsRoof.fimp==1
Vwater.VRoofSoilVeg(1,:,ittm)	=	repmat(0.*ParSoil.Roof.dz,1,1,1);
Owater.OwRoofSoilVeg(1,:,ittm)	=	0;
end

% Starting with dry soil if ground is fully impervious
if FractionsGround.fimp==1
Vwater.VGroundSoilImp(1,:,ittm)		=	repmat(0.*ParSoil.Ground.dz,1,1,1);
Owater.OwGroundSoilImp(1,:,ittm)	=	0;
Vwater.VGroundSoilBare(1,:,ittm)	=	repmat(0.*ParSoil.Ground.dz,1,1,1);
Owater.OwGroundSoilBare(1,:,ittm)	=	0;
Vwater.VGroundSoilVeg(1,:,ittm)		=	repmat(0.*ParSoil.Ground.dz,1,1,1);
Owater.OwGroundSoilVeg(1,:,ittm)	=	0;
end

% to save initial soil moisture for back-computation of water budget
OwaterInitial.OwRoofSoilVeg(1,:,ittm)     = Owater.OwRoofSoilVeg(1,:,ittm);
OwaterInitial.OwGroundSoilImp(1,:,ittm)   = Owater.OwGroundSoilImp(1,:,ittm);
OwaterInitial.OwGroundSoilBare(1,:,ittm)  = Owater.OwGroundSoilBare(1,:,ittm);
OwaterInitial.OwGroundSoilVeg(1,:,ittm)   = Owater.OwGroundSoilVeg(1,:,ittm);
OwaterInitial.OwGroundSoilTot(1,:,ittm)   = Owater.OwGroundSoilTot(1,:,ittm);


for ittn	= 1:n  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SunPosition,MeteoData,HumidityAtm,Anthropogenic,HVACSchedule,location,ParCalculation]...
	=feval(strcat('data_functions.UEHMForcingData_',Name_SiteFD),MeteoDataRaw,ittn,SoilPotW);


%% Initialize variables
[Gemeotry_m,ParTree,geometry,FractionsRoof,FractionsGround,...
	WallLayers,ParSoilRoof,ParSoilGround,ParInterceptionTree,...
	PropOpticalRoof,PropOpticalGround,PropOpticalWall,PropOpticalTree,...
	ParThermalRoof,ParThermalGround,ParThermalWall,ParThermalTree,...
	ParVegRoof,ParVegGround,ParVegTree,Person,...
    PropOpticalIndoors,ParThermalBulidingInt,ParWindows,ParHVAC,BEM_on]...
    =feval(strcat('data_functions.Data_UEHM_site_',Name_Site),MeteoData,ittm,LAI_TimeSeries);  

% Create variables with the values from the previous time step
[TempVec_ittm,Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,...
    CiCO2Leaf_ittm,Humidity_ittm,TempDamp_ittm,Runon_ittm,Qinlat_ittm,TempVecB_ittm,...
    ]=CreateVariablesOfPreviousTimeStep(ittm,ittn,TempVec,...
    Int,ExWater,Vwater,Owater,SoilPotW,CiCO2Leaf,Humidity,TempDamp,Runon,Qinlat,TempVecB,Results2m);

[TempVec_ittm2Ext,Humidity_ittm2Ext,TempVecB_ittm2Ext,Meteo_ittm]...
    =CreateVariablesOfPrevious2TimeStep(ittm,ittn,TempVec,Humidity,TempVecB,MeteoDataRaw);

if RESPreCalc == 1 || fconvPreCalc==1
    [fconv,rsRoofPreCalc,rsGroundPreCalc,rsTreePreCalc]=resistance_functions.PrecalculateForFasterNumericalSolution(ittn,ittm,...
             TempVec_ittm,Humidity_ittm,ParVegGround,SoilPotW_ittm,CiCO2Leaf_ittm,...
             MeteoData,HumidityAtm,geometry,FractionsGround,ParTree,...
		     PropOpticalGround,PropOpticalWall,PropOpticalTree,ParVegTree,...
             SunPosition,ViewFactor,ParWindows,BEM_on,ParVegRoof,...
             PropOpticalRoof,FractionsRoof,RES,Gemeotry_m);
else
    fconv = NaN; rsRoofPreCalc = NaN; rsGroundPreCalc = NaN; rsTreePreCalc = NaN;
end

% Turn on/off AC based on outdoor conditions
[ParHVAC,ParHVACorig]=BuildingEnergyModel.AC_HeatingTurnOnOff(ParHVAC,TempVecB_ittm,TempVec_ittm,Humidity_ittm,...
    MeteoData,Gemeotry_m,BEM_on);


for HVACittm = 1:2
% Only if itteration is 2    
if BEM_on==1 && HVACittm==2

    if ParHVACorig.ACon==0 && ParHVACorig.Heatingon==0
        break
    end

   if EnergyUse.EnergyForAC_H>-10^-6 && EnergyUse.EnergyForAC_LE>-10^-6 && EnergyUse.EnergyForHeating>-10^-6
        if ParHVAC.ACon==1 && round(TempVecB.Tbin(ittn,1,ittm),4)<(ParHVAC.TsetpointCooling+0.01) && round(TempVecB.qbin(ittn,1,ittm),8)<(ParHVAC.q_RHspCooling+10^-6)
            break
        elseif ParHVAC.Heatingon==1 && round(TempVecB.Tbin(ittn,1,ittm),4)>(ParHVAC.TsetpointHeating-0.01)
            break
        end
    end

    % Switch on or off HVAC based on temperature and humidity thresholds
    % (Humidity threshold is currently turned off as we assume that AC is
    % switched on mostly for cooling purposes and not for dehumidifaction
    if ParHVACorig.ACon==1 && round(TempVecB.Tbin(ittn,1,ittm),4)>(ParHVAC.TsetpointCooling+0.01) %|| round(TempVecB.qbin(ittn,1,ittm),8)>(ParHVAC.q_RHspCooling+10^-6)
        % Switch on AC based on exceedance of set-point temperature or humidity
        ParHVAC.ACon        = 1;
        ParHVAC.AC_onCool   = 1;
        ParHVAC.AC_onDehum  = 1;
        ParHVAC.Heatingon   = 0;
        ParHVAC.MasterOn    = 1;
        if round(TempVecB.qbin(ittn,1,ittm),8)<(ParHVAC.q_RHspCooling+10^-6) && round(EnergyUse.EnergyForAC_LE,1)==0
            ParHVAC.AC_onDehum = 0;
        elseif round(TempVecB.Tbin(ittn,1,ittm),4)<(ParHVAC.TsetpointCooling+0.01) && round(EnergyUse.EnergyForAC_H,1)==0
            ParHVAC.AC_onCool = 0;
        end
    elseif ParHVACorig.ACon==1 && EnergyUse.EnergyForAC_H>0 && round(TempVecB.qbin(ittn,1,ittm),8)>(ParHVAC.q_RHspCooling+10^-6)
        ParHVAC.ACon        = 1;
        ParHVAC.AC_onCool   = 1;
        ParHVAC.AC_onDehum  = 1;
        ParHVAC.MasterOn    = 1;
    elseif ParHVACorig.Heatingon==1 && round(TempVecB.Tbin(ittn,1,ittm),4)<(ParHVAC.TsetpointHeating-0.01) && round(EnergyUse.EnergyForHeating)==0
        % Switch on heating based on not reaching the set-point temperature
        ParHVAC.ACon        = 0;
        ParHVAC.AC_onCool   = 0;
        ParHVAC.AC_onDehum  = 0;
        ParHVAC.Heatingon   = 1;
        ParHVAC.MasterOn    = 1;
    end

    % Switch of HVAC because of negative energy consumption
    if ParHVACorig.ACon==1 && EnergyUse.EnergyForAC_H<-10^-6 || EnergyUse.EnergyForAC_LE<-10^-6
        if EnergyUse.EnergyForAC_H<-10^-6 && EnergyUse.EnergyForAC_LE<-10^-6
            ParHVAC.ACon        = 0;
            ParHVAC.AC_onCool   = 0;
            ParHVAC.AC_onDehum  = 0;
            ParHVAC.MasterOn    = 1;
        elseif EnergyUse.EnergyForAC_H<-10^-6 && EnergyUse.EnergyForAC_LE>10^-6
            ParHVAC.ACon        = 1;
            ParHVAC.AC_onCool   = 0;
            ParHVAC.AC_onDehum  = 1;
            ParHVAC.MasterOn    = 1;
        elseif EnergyUse.EnergyForAC_H>10^-6 && EnergyUse.EnergyForAC_LE<-10^-6
            ParHVAC.ACon        = 1;
            ParHVAC.AC_onCool   = 1;
            ParHVAC.AC_onDehum  = 0;
            ParHVAC.MasterOn    = 1;
        elseif EnergyUse.EnergyForAC_H<-10^-6 && EnergyUse.EnergyForAC_LE==0
            ParHVAC.ACon        = 0;
            ParHVAC.AC_onCool   = 0;
            ParHVAC.AC_onDehum  = 0;
            ParHVAC.MasterOn    = 1;
        elseif EnergyUse.EnergyForAC_LE<-10^-6 && EnergyUse.EnergyForAC_H==0
            ParHVAC.ACon        = 0;
            ParHVAC.AC_onCool   = 0;
            ParHVAC.AC_onDehum  = 0;
            ParHVAC.MasterOn    = 1;
        end
    elseif ParHVACorig.Heatingon==1 && EnergyUse.EnergyForHeating<-10^-6
            ParHVAC.Heatingon   = 0;
            ParHVAC.MasterOn    = 1;
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
% Solve energy budget equations for canyon, roof, and building interior
[Ttot,fval,exitflag]=fSolver_Tot(TempVec_ittm,TempVecB_ittm,Humidity_ittm,MeteoData,...
		Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,CiCO2Leaf_ittm,TempDamp_ittm,ViewFactor,...
		Gemeotry_m,ParTree,geometry,FractionsGround,FractionsRoof,...
		WallLayers,ParSoilGround,ParInterceptionTree,...
		PropOpticalGround,PropOpticalWall,PropOpticalTree,...
		ParThermalGround,ParThermalWall,ParVegGround,ParVegTree,...
        ParSoilRoof,PropOpticalRoof,ParThermalRoof,ParVegRoof,...
		SunPosition,HumidityAtm,Anthropogenic,ParCalculation,...
        PropOpticalIndoors,ParHVAC,ParThermalBulidingInt,ParWindows,BEM_on,...
        TempVec_ittm2Ext,Humidity_ittm2Ext,TempVecB_ittm2Ext,Meteo_ittm,...
        RESPreCalc,fconvPreCalc,fconv,rsRoofPreCalc,rsGroundPreCalc,rsTreePreCalc,HVACSchedule);


% Assign output temperatures and humidity
TempVec.TRoofImp(ittn,1,ittm)		=	Ttot(1,1);
TempVec.TRoofVeg(ittn,1,ittm)		=	Ttot(1,2);
TempVec.TRoofIntImp(ittn,1,ittm)	=	Ttot(1,3);
TempVec.TRoofIntVeg(ittn,1,ittm)	=	Ttot(1,4);

TempVec.TGroundImp(ittn,1,ittm)		=	Ttot(1,5);
TempVec.TGroundBare(ittn,1,ittm)	=	Ttot(1,6);
TempVec.TGroundVeg(ittn,1,ittm)		=	Ttot(1,7);
TempVec.TWallSun(ittn,1,ittm)		=	Ttot(1,8);
TempVec.TWallShade(ittn,1,ittm)		=	Ttot(1,9);
TempVec.TTree(ittn,1,ittm)			=	Ttot(1,10);
TempVec.TWallIntSun(ittn,1,ittm)	=	Ttot(1,11);
TempVec.TWallIntShade(ittn,1,ittm)	=	Ttot(1,12);
TempVec.TCanyon(ittn,1,ittm)		=	Ttot(1,13);
Humidity.CanyonSpecific(ittn,1,ittm)=	Ttot(1,14);

TempVecB.Tceiling(ittn,1,ittm)	    =	Ttot(1,15);
TempVecB.Tinwallsun(ittn,1,ittm)    =	Ttot(1,16);
TempVecB.Tinwallshd(ittn,1,ittm)    =	Ttot(1,17);
TempVecB.Twindows(ittn,1,ittm)     =	Ttot(1,18);
TempVecB.Tinground(ittn,1,ittm)     =	Ttot(1,19);
TempVecB.Tintmass(ittn,1,ittm)      =	Ttot(1,20);
TempVecB.Tbin(ittn,1,ittm)		    =	Ttot(1,21);
TempVecB.qbin(ittn,1,ittm)		    =	Ttot(1,22);

%--------------------------------------------------------------------------
TR(1,1) =   TempVec.TRoofImp(ittn,1,ittm);
TR(1,2) =   TempVec.TRoofVeg(ittn,1,ittm);
TR(1,3) =   TempVec.TRoofIntImp(ittn,1,ittm);
TR(1,4) =   TempVec.TRoofIntVeg(ittn,1,ittm);

TC(:,1) =	TempVec.TGroundImp(ittn,1,ittm);
TC(:,2) =	TempVec.TGroundBare(ittn,1,ittm);
TC(:,3) =	TempVec.TGroundVeg(ittn,1,ittm);
TC(:,4) =	TempVec.TWallSun(ittn,1,ittm);
TC(:,5) =	TempVec.TWallShade(ittn,1,ittm);
TC(:,6) =	TempVec.TTree(ittn,1,ittm);
TC(:,7) =	TempVec.TWallIntSun(ittn,1,ittm);
TC(:,8) =	TempVec.TWallIntShade(ittn,1,ittm);
TC(:,9) =	TempVec.TCanyon(ittn,1,ittm);
TC(:,10)=	Humidity.CanyonSpecific(ittn,1,ittm);

TB(:,1) =	TempVecB.Tceiling(ittn,1,ittm);
TB(:,2) =	TempVecB.Tinwallsun(ittn,1,ittm);
TB(:,3) =	TempVecB.Tinwallshd(ittn,1,ittm);
TB(:,4) =	TempVecB.Twindows(ittn,1,ittm);
TB(:,5) =	TempVecB.Tinground(ittn,1,ittm);
TB(:,6) =	TempVecB.Tintmass(ittn,1,ittm);
TB(:,7) =	TempVecB.Tbin(ittn,1,ittm);
TB(:,8) =	TempVecB.qbin(ittn,1,ittm);

% Calculate Energy and Water fluxes output for Roof
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
		=EB_WB_roof(TR,TB,TempVec_ittm,MeteoData,...
		Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,CiCO2Leaf_ittm,Runon_ittm,...
		Gemeotry_m,FractionsRoof,ParSoilRoof,PropOpticalRoof,ParThermalRoof,ParVegRoof,...
		HumidityAtm,Anthropogenic,ParCalculation,BEM_on,...
        RESPreCalc,rsRoofPreCalc);

	
% Calculate Energy and Water flux output for Canyon
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
		raCanyontoAtm,raCanyontoAtmOrig,rap_can,rap_Htree_In,rb_HGround,rb_LGround,...
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
		DEi,Eimp_2m,Ebare_soil_2m,Eveg_int_2m,Eveg_soil_2m,TEveg_2m,Ecan_2m,dS_H_air,dS_LE_air]...
		=EB_WB_canyon(TC,TB,TempVec_ittm,Humidity_ittm,MeteoData,...
		Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,...
		CiCO2Leaf_ittm,TempDamp_ittm,Runon_ittm,Qinlat_ittm,ViewFactor,...
		Gemeotry_m,ParTree,geometry,FractionsGround,...
		WallLayers,ParSoilGround,ParInterceptionTree,...
		PropOpticalGround,PropOpticalWall,PropOpticalTree,...
		ParThermalGround,ParThermalWall,ParVegGround,ParVegTree,...
		SunPosition,HumidityAtm,Anthropogenic,ParCalculation,...
        TempVecB_ittm,G2Roof,PropOpticalIndoors,ParHVAC,ParThermalBulidingInt,ParWindows,BEM_on,...
        RESPreCalc,fconvPreCalc,fconv,rsGroundPreCalc,rsTreePreCalc,HVACSchedule);


SWRinWsun = SWRabs_t.SWRabsWallSunTransmitted;
SWRinWshd = SWRabs_t.SWRabsWallShadeTransmitted;

% Calculate output variables of building energy model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[HbuildIntc,LEbuildIntc,GbuildIntc,SWRabsBc,LWRabsBc,TDampGroundBuild,WasteHeat,EnergyUse,HumidBuilding,ParACHeat_t,YBuildInt]=...
    BuildingEnergyModel.EBSolver_BuildingOUTPUT(TC,TB,TempVecB_ittm,TempVec_ittm,Humidity_ittm,MeteoData,...
    SWRinWsun,SWRinWshd,G2Roof,G2WallSun,G2WallShade,TempDamp_ittm,SWRabs_t,...
    Gemeotry_m,PropOpticalIndoors,ParHVAC,ParCalculation,ParThermalBulidingInt,ParWindows,BEM_on,HVACSchedule);

end

% Calculate Outdoor Thermal comfort outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Tmrt,BoleanInSun,SWRdir_Person,SWRdir_in_top,SWRdir_in_bottom,...
	SWRdir_in_east,SWRdir_in_south,SWRdir_in_west,SWRdir_in_north,...
	SWRdiff_Person,LWR_Person]=MRT.MeanRadiantTemperature(SWRout_t,LWRout_t,MeteoData,ViewFactorPoint,...
	ParTree,ParVegTree,geometry,Gemeotry_m,SunPosition,Person);
	
[u_ZPerson]=resistance_functions.WindProfile_PointOutput(Person.HeightWind,...
	Gemeotry_m,ParVegTree,ParTree,MeteoData,FractionsGround,ParVegGround);	

[UTCI_approx]=OTC.UTCI_approx(T2m-273.15,RH_T2m.*100,Tmrt,u_ZPerson);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate urban average of energy and water fluxes
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
% Numerical solver success ------------------------------------------------
Solver.Success(ittn,:,ittm)		=	exitflag;
Solver.ValuesEB(ittn,:,ittm)    =	fval;
Solver.Tsolver(ittn,:,ittm)		=	Ttot;
Solver.YfunctionOutput(ittn,:,ittm) =	[Yroof, Ycanyon, YBuildInt];

% Results 2m temperature & humidity ---------------------------------------
for i=1:length(Results2mNames)
	Results2m.(cell2mat(Results2mNames(i)))(ittn,1,ittm) = eval(cell2mat(Results2mNames(i)));
end

% Results 2m energy fluxes ------------------------------------------------
for i=1:length(Results2mEnergyFluxNames)
	Results2mEnergyFlux.(cell2mat(Results2mEnergyFluxNames(i)))(ittn,1,ittm) = eval(cell2mat(Results2mEnergyFluxNames(i)));
end

% Mean radiant temperature ------------------------------------------------
for i=1:length(MeanRadiantTemperatureNames)
	MeanRadiantTemperature.(cell2mat(MeanRadiantTemperatureNames(i)))(ittn,1,ittm)		=	eval(cell2mat(MeanRadiantTemperatureNames(i)));
end

% Universal thermal comfort index (UTCI) ----------------------------------
UTCI(ittn,1,ittm)	=	UTCI_approx;

% Humidity ----------------------------------------------------------------
for i=1:6
	Humidity.(cell2mat(HumidityNames(i)))(ittn,1,ittm)=	HumidityCan.(cell2mat(HumidityNames(i)));
end

for i=7:12
	Humidity.(cell2mat(HumidityNames(i)))(ittn,1,ittm)=	HumidityAtm.(cell2mat(HumidityNames(i)));
end

% Shortwave radiation outdoor ---------------------------------------------
for i=1:3
	SWRabs.(cell2mat(SWRabsNames(i)))(ittn,1,ittm)	=	eval(cell2mat(SWRabsNames(i)));
	SWRin.(cell2mat(SWRinNames(i)))(ittn,1,ittm)	=	eval(cell2mat(SWRinNames(i)));
	SWRout.(cell2mat(SWRoutNames(i)))(ittn,1,ittm)	=	eval(cell2mat(SWRoutNames(i)));
	SWREB.(cell2mat(SWREBNames(i)))(ittn,1,ittm)	=	eval(cell2mat(SWREBNames(i)));
end
for i=4:size(SWRinNames,1)
	SWRin.(cell2mat(SWRinNames(i)))(ittn,1,ittm)	=	SWRin_t.(cell2mat(SWRinNames(i)))(1,1);
	SWRout.(cell2mat(SWRoutNames(i)))(ittn,1,ittm)	=	SWRout_t.(cell2mat(SWRoutNames(i)))(1,1);
	SWREB.(cell2mat(SWREBNames(i)))(ittn,1,ittm)	=	SWREB_t.(cell2mat(SWREBNames(i)))(1,1);
end

for i=4:size(SWRabsNames,1)
	SWRabs.(cell2mat(SWRabsNames(i)))(ittn,1,ittm)	=	SWRabs_t.(cell2mat(SWRabsNames(i)))(1,1);
end

% Longwave radiation outdoor ----------------------------------------------
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

% Sensible heat outdoor ---------------------------------------------------
for i=1:size(HfluxNames,1)
	Hflux.(cell2mat(HfluxNames(i)))(ittn,1,ittm)	=	eval(cell2mat(HfluxNames(i)));
end

% Latent heat outdoor -----------------------------------------------------
for i=1:size(LEfluxNames,1)
	LEflux.(cell2mat(LEfluxNames(i)))(ittn,1,ittm)	=	eval(cell2mat(LEfluxNames(i)));
end

% Ground heat flux outdoor ------------------------------------------------
for i=1:size(GfluxNames,1)
	Gflux.(cell2mat(GfluxNames(i)))(ittn,1,ittm)	=	eval(cell2mat(GfluxNames(i)));
end

% Heat Storage ------------------------------------------------------------
for i=1:size(dStorageNames,1)
	dStorage.(cell2mat(dStorageNames(i)))(ittn,1,ittm)	=	eval(cell2mat(dStorageNames(i)));
end

% Dampening temperature ---------------------------------------------------
for i=1:size(TempDampNames,1)
	TempDamp.(cell2mat(TempDampNames(i)))(ittn,1,ittm)	=	eval(cell2mat(TempDampNames(i)));
end

% Resistances for H and LE ------------------------------------------------
for i=1:size(RESNames,1)
	RES.(cell2mat(RESNames(i)))(ittn,1,ittm)	=	eval(cell2mat(RESNames(i)));
end

% Evapotranspiration flux -------------------------------------------------
for i=1:size(EfluxNames,1)
	Eflux.(cell2mat(EfluxNames(i)))(ittn,1,ittm)	=	eval(cell2mat(EfluxNames(i)));
end

% Runoff ------------------------------------------------------------------
for i=1:size(RunoffNames,1)
	Runoff.(cell2mat(RunoffNames(i)))(ittn,1,ittm)	=	eval(cell2mat(RunoffNames(i)));
end

% Runon -------------------------------------------------------------------
for i=1:size(RunonNames,1)
	Runon.(cell2mat(RunonNames(i)))(ittn,1,ittm)	=	eval(cell2mat(RunonNames(i)));
end

% Leakage at the bottom of the soil ---------------------------------------
for i=1:size(LeakageNames,1)
	Leakage.(cell2mat(LeakageNames(i)))(ittn,1,ittm)	=	eval(cell2mat(LeakageNames(i)));
end

% Interception and ponding water ------------------------------------------
for i=1:size(IntNames,1)
	Int.(cell2mat(IntNames(i)))(ittn,1,ittm)	=	eval(cell2mat(IntNames(i)));
end

% Change in interception and ponding water --------------------------------
for i=1:size(dInt_dtNames,1)
	dInt_dt.(cell2mat(dInt_dtNames(i)))(ittn,1,ittm)	=	eval(cell2mat(dInt_dtNames(i)));
end

% Infiltration ------------------------------------------------------------
for i=1:size(InfiltrationNames,1)
	Infiltration.(cell2mat(InfiltrationNames(i)))(ittn,1,ittm)=	eval(cell2mat(InfiltrationNames(i)));
end

% Water volumne in soil ---------------------------------------------------
VwaterNames	=	fieldnames(Vwater);
for i=1:size(VwaterNames,1)
	Vwater.(cell2mat(VwaterNames(i)))(ittn,:,ittm)=	eval(cell2mat(VwaterNames(i)));
end

% Change in water volumne in soil -----------------------------------------
dVwater_dtNames	=	fieldnames(dVwater_dt);
for i=1:size(dVwater_dtNames,1)
	dVwater_dt.(cell2mat(dVwater_dtNames(i)))(ittn,:,ittm)=	eval(cell2mat(dVwater_dtNames(i)));
end

% Water content in soil ---------------------------------------------------
OwaterNames	=	fieldnames(Owater);
for i=1:size(OwaterNames,1)
	Owater.(cell2mat(OwaterNames(i)))(ittn,:,ittm)=	eval(cell2mat(OwaterNames(i)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Water content in soil ---------------------------------------------------
OSwaterNames	=	fieldnames(OSwater);
for i=1:size(OSwaterNames,1)
	OSwater.(cell2mat(OSwaterNames(i)))(ittn,:,ittm)=	eval(cell2mat(OSwaterNames(i)));
end

% Lateral soil water flux -------------------------------------------------
for i=1:size(QinlatNames,1)
	Qinlat.(cell2mat(QinlatNames(i)))(ittn,:,ittm)=	eval(cell2mat(QinlatNames(i)));
end

% Extractable water from soil ---------------------------------------------
for i=1:size(ExWaterNames,1)
	ExWater.(cell2mat(ExWaterNames(i)))(ittn,:,ittm)=	eval(cell2mat(ExWaterNames(i)));
end

% Soil water potential ----------------------------------------------------
for i=1:size(SoilPotWNames,1)
	SoilPotW.(cell2mat(SoilPotWNames(i)))(ittn,1,ittm)=	eval(cell2mat(SoilPotWNames(i)));
end

% CO2 concentration within leaf -------------------------------------------
for i=1:size(CiCO2LeafNames,1)
	CiCO2Leaf.(cell2mat(CiCO2LeafNames(i)))(ittn,1,ittm)=	eval(cell2mat(CiCO2LeafNames(i)));
end

% Energy budget closure ---------------------------------------------------
for i=1:size(EBNames,1)
	EB.(cell2mat(EBNames(i)))(ittn,1,ittm)=	eval(cell2mat(EBNames(i)));
end

% Water budget closure roof -----------------------------------------------
WBRoofVegSoil	=	sum(WBRoofVegSoil);
for i=1:size(WBRoofNames,1)
	WBRoof.(cell2mat(WBRoofNames(i)))(ittn,:,ittm)=	eval(cell2mat(WBRoofNames(i)));
end

% Water budget closure individual -----------------------------------------
WBCanyonIndv.WB_In_tree(ittn,:,ittm)	=	WBIndv.WB_In_tree;
WBCanyonIndv.WB_In_gveg(ittn,:,ittm)	=	WBIndv.WB_In_gveg;
WBCanyonIndv.WB_In_gimp(ittn,:,ittm)	=	WBIndv.WB_In_gimp;
WBCanyonIndv.WB_In_gbare(ittn,:,ittm)	=	WBIndv.WB_In_gbare;
WBCanyonIndv.WB_Pond_gveg(ittn,:,ittm)	=	WBIndv.WB_Pond_gveg;
WBCanyonIndv.WB_Soil_gimp(ittn,:,ittm)	=	nansum(WBIndv.WB_Soil_gimp);
WBCanyonIndv.WB_Soil_gbare(ittn,:,ittm)	=	sum(WBIndv.WB_Soil_gbare);
WBCanyonIndv.WB_Soil_gveg(ittn,:,ittm)	=	sum(WBIndv.WB_Soil_gveg);

% Water budget closure canyon ---------------------------------------------
WBCanyonTotNames	=	fieldnames(WBCanyonTot);
for i=1:size(WBCanyonTotNames,1)
	WBCanyonTot.(cell2mat(WBCanyonTotNames(i)))(ittn,:,ittm)=	WBTot.(cell2mat(WBCanyonTotNames(i)));
end

% Wind speed outputs at different heights ---------------------------------
for i=1:size(WindNames,1)
	Wind.(cell2mat(WindNames(i)))(ittn,1,ittm)=	eval(cell2mat(WindNames(i)));
end

%--------------------------------------------------------------------------
% Building energy model
% Sensible heat flux building interior ------------------------------------
for i=1:length(HbuildIntNames)
	HbuildInt.(cell2mat(HbuildIntNames(i)))(ittn,1,ittm)=	HbuildIntc.(cell2mat(HbuildIntNames(i)));
end

% Latent heat flux building interior --------------------------------------
for i=1:length(LEbuildIntNames)
	LEbuildInt.(cell2mat(LEbuildIntNames(i)))(ittn,1,ittm)=	LEbuildIntc.(cell2mat(LEbuildIntNames(i)));
end

% Conductive heat flux building interior ----------------------------------
for i=1:length(GbuildIntNames)
	GbuildInt.(cell2mat(GbuildIntNames(i)))(ittn,1,ittm)=	GbuildIntc.(cell2mat(GbuildIntNames(i)));
end

% Absorbed shortwave radiation building interior --------------------------
for i=1:length(SWRabsBNames)
	SWRabsB.(cell2mat(SWRabsBNames(i)))(ittn,1,ittm)=	SWRabsBc.(cell2mat(SWRabsBNames(i)));
end

% Absorbed longwave radiation building interior ---------------------------
for i=1:length(LWRabsBNames)
	LWRabsB.(cell2mat(LWRabsBNames(i)))(ittn,1,ittm)=	LWRabsBc.(cell2mat(LWRabsBNames(i)));
end

% Waste heat emitted into canyon air --------------------------------------
for i=1:length(BEMWasteHeatNames)
	BEMWasteHeat.(cell2mat(BEMWasteHeatNames(i)))(ittn,:,ittm)=	WasteHeat.(cell2mat(BEMWasteHeatNames(i)));
end

% Humidity within the canyon --------------------------------------
for i=1:length(HumidityBuildingNames)
	HumidityBuilding.(cell2mat(HumidityBuildingNames(i)))(ittn,:,ittm)=	HumidBuilding.(cell2mat(HumidityBuildingNames(i)));
end

% Building energy usage ---------------------------------------------------
for i=1:length(BEMEnergyUseNames)
	BEMEnergyUse.(cell2mat(BEMEnergyUseNames(i)))(ittn,:,ittm)=	EnergyUse.(cell2mat(BEMEnergyUseNames(i)));
end

% BEMEnergyUse.EnergyForAC(ittn,:,ittm) = EnergyUse.EnergyForAC;
% BEMEnergyUse.EnergyForHeating(ittn,:,ittm) = EnergyUse.EnergyForHeating;

% LAI timeseries ----------------------------------------------------------
LAI_ts.LAI_R(ittn,1,ittm)	=	ParVegRoof.LAI;
LAI_ts.LAI_G(ittn,1,ittm)	=	ParVegGround.LAI;
LAI_ts.LAI_T(ittn,1,ittm)	=	ParVegTree.LAI;

% Anthropogenic heat and water input --------------------------------------
AnthropoNames = fieldnames(Anthropogenic);
for i=1:length(AnthropoNames)
	Anthropo.(cell2mat(AnthropoNames(i)))(ittn,:,ittm)=	Anthropogenic.(cell2mat(AnthropoNames(i)));
end

% Albedo ------------------------------------------------------------------
albedo_urban	=	geometry.wcanyon_norm*albedo_canyon + geometry.wroof_norm*PropOpticalRoof.albedo;

AlbedoOutput.TotalUrban(ittn,:,ittm)	=	albedo_urban;
AlbedoOutput.TotalCanyon(ittn,:,ittm)	=	albedo_canyon;
AlbedoOutput.Roof(ittn,:,ittm)			=	PropOpticalRoof.albedo;

% AC and heating switches time varying --------------------------------------
ParACHeatNames = fieldnames(ParACHeat_ts);
for i=1:length(ParACHeatNames)
	ParACHeat_ts.(cell2mat(ParACHeatNames(i)))(ittn,:,ittm)=	ParACHeat_t.(cell2mat(ParACHeatNames(i)));
end


if mod(ittn,10)==0
disp(strcat('iter=',num2str(ittn),' iterm=',num2str(ittm)));
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

ParSoilGround.Osat  =   ParSoil.Ground.Osat;
ParSoilGround.Ohy   =   ParSoil.Ground.Ohy;
ParSoilGround.O33   =   ParSoil.Ground.O33;
ParSoilGround.dz    =   ParSoil.Ground.dz;
ParSoilRoof.Osat    =   ParSoil.Roof.Osat;
ParSoilRoof.Ohy     =   ParSoil.Roof.Ohy;
ParSoilRoof.O33     =   ParSoil.Roof.O33;
ParSoilRoof.dz      =   ParSoil.Roof.dz;
%--------------------------------------------------------------------------
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

PropOpticalIndoors_Out(ittm)    =	PropOpticalIndoors;
ParThermalBulidingInt_Out(ittm) =	ParThermalBulidingInt;
ParWindows_Out(ittm)			=	ParWindows;
ParHVAC_Out(ittm)		        =	ParHVAC;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
ParHVAC_Out.Acon = ParHVACorig.ACon;
ParHVAC_Out.Heatingon = ParHVACorig.Heatingon;

SwitchOnFigure = 1;

ParHVAC= ParHVAC_Out;

% Plot and calculate radiation and energy balance
[WaterFluxRoof,WaterFluxCan,WaterFluxBuild,WaterFluxUrban]=WaterBalanceComponents(MeteoDataRaw,...
    Runon,Leakage,LEflux,dVwater_dt,OwaterInitial,Owater,dInt_dt,Int,Anthropo,...
    LEbuildInt,BEMWasteHeat,ParHVAC,...
    ParSoilRoof_Out,ParSoilGround_Out,ParCalculation_Out,FractionsRoof_Out,...
    FractionsGround_Out,geometry_Out,Gemeotry_m_Out,SwitchOnFigure,ittm);

UrbanClimateVariables(TempVec,UTCI,Results2m,MeteoDataRaw,MeanRadiantTemperature,...
    TempVecB,HumidityBuilding,...
    FractionsGround_Out,FractionsRoof_Out,ParTree_Out,SwitchOnFigure,ittm,BEM_on);

[EnergyFluxUrban,EnergyFluxCan,EnergyFluxRoof]=PlanAreaEnergyBalanceCalculation(ViewFactor,MeteoDataRaw,...
    SWRin,SWRout,SWRabs,LWRin,LWRout,LWRabs,LEflux,Hflux,Gflux,BEMWasteHeat,...
    dStorage,GbuildInt,HbuildInt,LEbuildInt,SWRabsB,LWRabsB,...
    geometry_Out,FractionsGround_Out,PropOpticalRoof_Out,Anthropo,Gemeotry_m_Out,SwitchOnFigure,ittm,BEM_on);


if OutputsToSave==1
% Option 1: Essential energy flux and climate outputs
%--------------------------------------------------------------------------
save(['Calculation',NameOutput],'NameOutput','Solver','n','m','Name_Site','MeteoDataRaw',...
    'TempVec','Humidity','Results2m','TempVecB','HumidityBuilding',...
    'MeanRadiantTemperature','UTCI','Wind',...
    'EnergyFluxUrban','EnergyFluxCan','EnergyFluxRoof',...
    'WaterFluxRoof','WaterFluxBuild','WaterFluxCan','WaterFluxUrban',...
    'BEMEnergyUse','BEMWasteHeat',...
    ...
	'Gemeotry_m_Out','ParTree_Out','geometry_Out','FractionsRoof_Out','FractionsGround_Out',...
	'WallLayers_Out','ParSoilRoof_Out','ParSoilGround_Out',...
	'PropOpticalRoof_Out','PropOpticalGround_Out','PropOpticalWall_Out','PropOpticalTree_Out',...
	'ParThermalRoof_Out','ParThermalGround_Out','ParThermalWall_Out','ParThermalTree_Out',...
    'ParCalculation_Out','ParVegRoof_Out','ParVegGround_Out','ParVegTree_Out',...
    'PropOpticalIndoors_Out','ParThermalBulidingInt_Out','ParWindows_Out','ParHVAC_Out',...
    'Zatm','AlbedoOutput','ViewFactor','ParACHeat_ts','BEM_on');

elseif OutputsToSave==2
% Option 2: Extended energy flux and climate outputs
%--------------------------------------------------------------------------
save(['Calculation',NameOutput],'NameOutput','Solver','n','m','Name_Site','MeteoDataRaw',...
    'TempVec','Humidity','Results2m','TempVecB','HumidityBuilding',...
    'MeanRadiantTemperature','UTCI','Wind',...
    'EnergyFluxUrban','EnergyFluxCan','EnergyFluxRoof',...
    'WaterFluxRoof','WaterFluxBuild','WaterFluxCan','WaterFluxUrban',...
    ....
    'BEMEnergyUse','BEMWasteHeat',...
    ...
    'SWRabs','LWRabs','Hflux','LEflux','Gflux','dStorage',...
    'Eflux','Runoff','Runon','Leakage','Int','dInt_dt','Infiltration',...
    'Vwater','dVwater_dt','Owater','SoilPotW','LAI_ts',...
	'TempDamp','Anthropo',...
    'HbuildInt','LEbuildInt','GbuildInt','SWRabsB','LWRabsB',...
    ...
	'Gemeotry_m_Out','ParTree_Out','geometry_Out','FractionsRoof_Out','FractionsGround_Out',...
	'WallLayers_Out','ParSoilRoof_Out','ParSoilGround_Out',...
	'PropOpticalRoof_Out','PropOpticalGround_Out','PropOpticalWall_Out','PropOpticalTree_Out',...
	'ParThermalRoof_Out','ParThermalGround_Out','ParThermalWall_Out','ParThermalTree_Out',...
    'ParCalculation_Out','ParVegRoof_Out','ParVegGround_Out','ParVegTree_Out',...
    'PropOpticalIndoors_Out','ParThermalBulidingInt_Out','ParWindows_Out','ParHVAC_Out',...
    'Zatm','AlbedoOutput','ViewFactor','ParACHeat_ts','BEM_on');

else
% Option 3: Extended outputs 
%--------------------------------------------------------------------------
save(['Calculation',NameOutput],'NameOutput','Solver','TempVec','Humidity','SWRabs','SWRin','SWRout','SWREB','LWRabs','LWRin','LWRout',...
	'LWREB','Hflux','LEflux','Gflux','dStorage','RES','Eflux','Runoff','Runon','Leakage',...
	'Int','dInt_dt','Infiltration','Vwater','dVwater_dt','Owater',...
    'OwaterInitial','OSwater','ExWater','SoilPotW',...
	'CiCO2Leaf','WBRoof','WBCanyonIndv','WBCanyonTot','EB','Wind','TempDamp','Qinlat','Results2m',...
	'n','m','Name_Site','MeteoDataRaw','Anthropo',...
	'Gemeotry_m_Out','ParTree_Out','geometry_Out','FractionsRoof_Out','FractionsGround_Out',...
	'WallLayers_Out','ParSoilRoof_Out','ParSoilGround_Out','ParInterceptionTree_Out',...
	'PropOpticalRoof_Out','PropOpticalGround_Out','PropOpticalWall_Out','PropOpticalTree_Out',...
	'ParThermalRoof_Out','ParThermalGround_Out','ParThermalWall_Out','ParThermalTree_Out','ParCalculation_Out',...
	'ParVegRoof_Out','ParVegGround_Out','ParVegTree_Out',...
    'PropOpticalIndoors_Out','ParThermalBulidingInt_Out','ParWindows_Out','ParHVAC_Out',...
    'LAI_ts','Results2mEnergyFlux','MeanRadiantTemperature','Zatm','UTCI',...
    'AlbedoOutput','ViewFactor','EnergyFluxUrban','EnergyFluxCan','EnergyFluxRoof',...
    'WaterFluxRoof','WaterFluxBuild','WaterFluxCan','WaterFluxUrban',...
    'TempVecB','HbuildInt','LEbuildInt','GbuildInt','SWRabsB','LWRabsB',...
    'BEMEnergyUse','BEMWasteHeat','HumidityBuilding','ParACHeat_ts','BEM_on','ittm');

end
