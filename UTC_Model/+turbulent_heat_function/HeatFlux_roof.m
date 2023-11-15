function[Hroof_imp,Hroof_veg,Eroof_imp,Eroof_veg,Eroof_ground,Eroof_soil,TEroof_veg,...
	LEroof_imp,LEroof_veg,LEroof_ground,LEroof_soil,LTEroof_veg,...
	Ci_sun_roof,Ci_shd_roof,ra,rb,rap_L,r_soil,rs_sun,rs_shd]...
	=HeatFlux_roof(TemperatureR,TempVec_ittm,MeteoData,HumidityAtm,ParVegRoof,FractionsRoof,Gemeotry_m,...
	ParSoilRoof,ParCalculation,SoilPotW_ittm,Owater_ittm,Vwater_ittm,ExWater_ittm,Int_ittm,CiCO2Leaf_ittm,...
	SWRabs_dir,SWRabs_diff,RESPreCalc,rsRoofPreCalc)


% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hroof_imp			=	sensible heat from impervious area [W/m^2]
% Hroof_veg			=	sensible heat from vegetated area [W/m^2]
% Eroof_imp			=	Evaporation from ponding water on impervious area [kg/m^2 s]
% Eroof_veg			=	Evaporation from intercepted water on vegetation [kg/m^2 s]
% Eroof_ground		=	Evaporation from ponding water on ground under vegetation [kg/m^2 s]
% Eroof_soil		=	Evaporation from first soil layer [kg/m^2 s]
% TEroof_veg		=	Transpiration from vegetation transpiration [kg/m^2 s]
% LEroof_imp		=	Latent heat from ponding water on impervious area [W/m^2]
% LEroof_veg		=	Latent heat from intercepted water on vegetation [W/m^2]
% LEroof_ground		=	Latent heat from ponding water on ground under vegetation [W/m^2]
% LEroof_soil		=	Latent heat from first soil layer [W/m^2]
% LTEroof_veg		=	Latent heat from vegetation transpiration [W/m^2]
% Ci_sun_roof		=	Leaf Interior  CO2 concentration [umolCO2/mol]
% Ci_shd_roof		=	Leaf Interior  CO2 concentration [umolCO2/mol]
% ra				=	Aerodynamic resistance [s/m] 
% rb_L				=	Leaf boundary resistance [s/m] 
% rap_L				=	Undercanopy resistance [s/m] 
% r_soil			=	Soil resistance [s/m] 
% rs_sun			=	Stomata resistance of sunlit vegetation [s/m] 
% rs_shd			=	Stomata resistance of shaded vegetation [s/m] 

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Troof_imp			=	Temperature of the roof impervious [K]
% Troof_veg			=	Temperature of the roof vegetated [K]
% H					=	Canyon (building) height (m)
% fveg_roof			=	Partitioning roof vegetation surface(-)
% fimp_roof			=	Partitioning roof impervious(-)
% LAI_roof			=	Leaf area index for the roof vegetation (-)
% SAI_roof			=	Stem area index for the roof vegetation (-)
% hc_roof			=	canopy height roof vegetation
% d_leaf_roof		=	Leaf dimension of roof vegetation [cm]
% Zatm				=	Atmospheric reference height [m]
% Tatm				=	Air Temperature at atmospheric reference level [K]
% Uatm				=	Wind speed at atmospheric reference level [m/s]
% Pre				=	air pressure [Pa]. Carefull, used to be [hPa - mbar] in T&C
% Catm_CO2			=	[ppm]-[umolCO2/mol] Atmospheric CO2 concentration 2017
% esat_Tatm			=	vapor pressure saturation at Tatm [Pa]
% ea				=	vapor pressure [Pa]
% SWRabs_dir		=	direct incoming shortwave radiation [W/m^2]
% SWRabs_diff		=	diffuse incoming shortwave radiation [W/m^2]
% Catm_O2			=	Intercellular Partial Pressure Oxygen [umolO2/mol]
% FI_roof			=	Intrinsec quantum Efficiency [umolCO2/umolPhotons]
% Do_roof			=	[Pa] Empirical coefficient for the role of vapor pressure in the biochemical model of photosynthesis
% a1_roof			=	[-] Empirical parameter connecting stomatal aperture and net assimilaton
% go_roof			=	[mol / s m^2] minimal Stomatal Conductance
% CT_roof			=	'CT' == 3  'CT' ==  4 Photosyntesis Typology for Plants, Photosynthetic pathway C3 or C4
% DSE_roof			=	[kJ/mol] Activation Energy - Plant Dependent, Activation Energy in Photosynthesis for Rubisco Capacity
% Ha_roof			=	[kJ / mol K]  entropy factor - Plant Dependent, Activation energy.
% gmes_roof			=	[mol CO2/s m2] Mesophyll conductance, not used
% rjv_roof			=	[?mol Eq/ ?mol CO2] Scaling factor between Jmax and Vmax
% Kopt_roof			=	[-] optical depth of direct beam perunit plant area ???
% Knit_roof			=	[-] Canopy nitrogen decay coefficient
% Vmax_roof			=	[?mol CO2/ m2 s] Maximum Rubisco capacity at 25°C leaf level
% mSl_roof			=	
% e_rel_roof		=	[-] Relative Efficiency of the photosynthesis apparatus due to Age/Day-length
% e_relN_roof		=	[-] Relative efficiency of the photosynthesis apparatus due to N limitations
% Psi_sto_00_roof	=	[MPa]  Water Potential at PLCs loss conductivity
% Psi_sto_50_roof	=	[MPa]  Water Potential at 50% loss conductivity
% Sl_roof			=	0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
% Sp_In_roof		=	specific water retained by a vegetated surface [mm m^2 VEG area m^-2 plant area]
% Rrootl			=	Root length index [m root / m^2 PFT]
% PsiL50			=	[MPa]  Water Potential at 50% loss conductivity
% PsiX50			=	Water potential at 50 of xylem hyd
% Pcla				=	Fraction of clay in the soil [-]
% Psan				=	Fraction of sand in the soil [-]
% Porg				=	Fraction of organic material in the soil [-]
% Kfc				=	Conductivity at field capacity [mm/h]
% Phy				=	Suction at the residual/hygroscopic water content [kPa]
% SPAR				=	SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
% Kbot				=	[mm/h] Conductivity at the bedrock layer
% CASE_ROOT			=	Type of Root Profile
% ZR95				=	Root depth 95 percentile [mm]
% ZR50				=	Root depth 50 percentile [mm]
% ZRmax				=	Maximum Root depth [mm]
% Zs				=	soil layer discretization [mm]
% In_imp_tm1		=	Interception on impervious area at the previous time step [mm]
% In_veg_tm1		=	Interception on plant surfaces at the previous time step [mm]
% In_ground_tm1		=	Interception on ground at the previous time step [mm]
% Vtm1				=	Water volume in the soil layers at prevous time step [mm]
% Exwat_tm1			=	Max extractable water for plants in soil
% Ci_sun_tm1		=	Leaf Interior  CO2 concentration [umolCO2/mol] at prevous time step
% Ci_shd_tm1		=	Leaf Interior  CO2 concentration [umolCO2/mol] at prevous time step
% Psi_ltm1			=	Soil water potential for plants in soil
% Otm1				=	Water content in the soil layers at previous time step [-]
% dth				=	time step of calculation [h]
% row				=	density of water [kg/m^3]

Troof_imp		=	TemperatureR(1,1);
Troof_veg		=	TemperatureR(1,2);
Troof_veg_tm1   =	TempVec_ittm.TRoofVeg;
Tatm			=	MeteoData.Tatm;
Pre				=	MeteoData.Pre;
ea				=	MeteoData.ea;
Zatm			=	MeteoData.Zatm;
Uatm			=	MeteoData.Uatm;
Catm_O2			=	MeteoData.Catm_O2;
Catm_CO2		=	MeteoData.Catm_CO2;
esat_Tatm		=	HumidityAtm.AtmVapourPreSat;
H				=	Gemeotry_m.Height_canyon;
fveg_roof		=	FractionsRoof.fveg;
fimp_roof		=	FractionsRoof.fimp;
LAI_roof		=	ParVegRoof.LAI;
SAI_roof		=	ParVegRoof.SAI;
hc_roof			=	ParVegRoof.hc;
d_leaf_roof		=	ParVegRoof.d_leaf;
Kopt_roof		=	ParVegRoof.Kopt;
Knit_roof		=	ParVegRoof.Knit;
Psi_sto_50_roof	=	ParVegRoof.Psi_sto_50;
Psi_sto_00_roof	=	ParVegRoof.Psi_sto_00;
CT_roof			=	ParVegRoof.CT;
Vmax_roof		=	ParVegRoof.Vmax;
DSE_roof		=	ParVegRoof.DSE;
Ha_roof			=	ParVegRoof.Ha;
FI_roof			=	ParVegRoof.FI;
Do_roof			=	ParVegRoof.Do;
a1_roof			=	ParVegRoof.a1;
go_roof			=	ParVegRoof.go;
e_rel_roof		=	ParVegRoof.e_rel;
e_relN_roof		=	ParVegRoof.e_relN;
gmes_roof		=	ParVegRoof.gmes;
rjv_roof		=	ParVegRoof.rjv;
mSl_roof		=	ParVegRoof.mSl;
Sl_roof			=	ParVegRoof.Sl;
CASE_ROOT		=	ParVegRoof.CASE_ROOT;
ZR95			=	ParVegRoof.ZR95;
ZR50			=	ParVegRoof.ZR50;
ZRmax			=	ParVegRoof.ZRmax;
Pcla			=	ParSoilRoof.Pcla;
Psan			=	ParSoilRoof.Psan;
Porg			=	ParSoilRoof.Porg;
Kfc				=	ParSoilRoof.Kfc;
Phy				=	ParSoilRoof.Phy;
SPAR			=	ParSoilRoof.SPAR;
Kbot			=	ParSoilRoof.Kbot;
Sp_In_roof		=	ParSoilRoof.Sp_In;
Zs				=	ParSoilRoof.Zs;
dth				=	ParCalculation.dth;
row				=	ParCalculation.row;

Psi_ltm1		=	SoilPotW_ittm.SoilPotWRoofVeg_L;
Otm1			=	Owater_ittm.OwRoofSoilVeg;
Vtm1			=	Vwater_ittm.VRoofSoilVeg;
Exwat_tm1		=	ExWater_ittm.ExWaterRoofVeg_L;
In_imp_tm1		=	Int_ittm.IntRoofImp;
In_veg_tm1		=	Int_ittm.IntRoofVegPlant;
In_ground_tm1	=	Int_ittm.IntRoofVegGround;
Ci_sun_tm1		=	CiCO2Leaf_ittm.CiCO2LeafRoofVegSun;
Ci_shd_tm1		=	CiCO2Leaf_ittm.CiCO2LeafRoofVegShd;


%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cp_atm		=	1005+(((Tatm-273.15)+23.15)^2)/3364;		% specific heat air  [J/kg K]
rho_atm		=	(Pre/(287.04*Tatm))*(1-(ea/Pre)*(1-0.622));	% dry air density at atmosphere [kg/m^3]
q_atm		=	0.622*ea/(Pre-0.378*ea);					% Specifc humidity of air at reference height []
L_heat		=	1000*(2501.3 - 2.361*(Tatm-273.15));				% Latent heat vaporization/condensaition [J/kg]

% Vapor pressure saturation and specific humidity at esat
esat_T_rimp	=	611*exp(17.27*(Troof_imp-273.16)/(237.3+(Troof_imp-273.16)));	% vapor pressure saturation at Troof_imp [Pa]
qsat_T_rimp	=	(0.622*esat_T_rimp)/(Pre-0.378*esat_T_rimp);					% Saturated specific humidity at Troof_imp []
esat_T_rveg	=	611*exp(17.27*(Troof_veg-273.16)/(237.3+(Troof_veg-273.16)));	% vapor pressure saturation at Troof_veg [Pa]
qsat_T_rveg	=	(0.622*esat_T_rveg)/(Pre-0.378*esat_T_rveg);					% Saturated specific humidity at Troof_veg []
Troof		=	fveg_roof*Troof_veg+fimp_roof*Troof_imp;						% Average temperature of the roof [K]
esat_roof	=	fveg_roof*esat_T_rveg+fimp_roof*esat_T_rimp;						% Average vapour pressure of the roof [K]

% Parameters for stomata resistance
Citm1_sun	=	Ci_sun_tm1;		% Ci = Leaf Interior  CO2 concentration [umolCO2/mol]
Citm1_shd	=	Ci_shd_tm1;		% Ci = Leaf Interior  CO2 concentration [umolCO2/mol]
Ds_atm		=	esat_Tatm-ea;	% Ds = Vapor Pressure Deficit [Pa]
Oa			=	Catm_O2;		% Oa Intercellular Partial Pressure Oxygen [umolO2/mol]

% Partitioning of radiation into sunlit and shaded area
Fsun			=	(1.0 - exp(-Kopt_roof*(LAI_roof)))/(Kopt_roof*(LAI_roof));
Fsun(Fsun<0.01)	=	0; 
Fsun(Fsun>1)	=	1;
Fshd			=	1- Fsun;
PAR_sun			=	SWRabs_dir+Fsun*SWRabs_diff;	% [W/m^2] absorbed direct and diffuse shortwave radiation of the sunlit surface
PAR_shd			=	Fshd*SWRabs_diff;				% [W/m^2] absorbed direct and diffuse shortwave radiation of the shaded surface 

% Fraction of vegetation covered by intercepted water (Deardorff(1978))
In_max_veg		=	Sp_In_roof*(LAI_roof+SAI_roof);
dw_veg_roof		=	min(1,(In_veg_tm1/In_max_veg)^(2/3));        % Fraction of vegetation covered by intercepted water

%% Calculation of resistances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fveg_roof==0
    hc_roof=0;
    LAI_roof=0;
    d_leaf_roof=0;
end

if fimp_roof==0
    Croof=0;
else
	Croof=1;
end

% Roughness height and length
[zom,zoh,~,~,disp_h,~,~,~,~,~,~,~]=resistance_functions.Urban_roughness(hc_roof,0,0,0,Croof);

% Wind profile
[u_ref_und,u_Hveg]=resistance_functions.WindProfile_Roof(H,Zatm,Uatm,zom,disp_h,hc_roof,2);

[ra]=resistance_functions.Aerodynamic_Resistence(Tatm-273.15,Troof-273.15,Pre/100,Zatm-H,disp_h,zom,zoh,Uatm,ea,esat_roof);
%[ra]=Aerodynamic_Resistence(Ta,Ts,Pre,zatm,disp_h,zom,zoh,Ws,ea,es)

[rb]=resistance_functions.Leaf_Boundary_Resistence(u_Hveg,Troof_veg-273.15,Tatm-273,hc_roof,d_leaf_roof,LAI_roof,Zatm,disp_h,zom);
% [rb]=Leaf_Boundary_Resistence(Ws,Ts,Ta,hc,d_leaf,LAI,zatm,disp_h,zom)


if RESPreCalc==1
    rs_sun=rsRoofPreCalc.rs_sun;  
    rs_shd=rsRoofPreCalc.rs_shd; 
    Ci_sun=rsRoofPreCalc.Ci_sun; 
    Ci_shd=rsRoofPreCalc.Ci_shd;
else
    % Stomatal resistance,roughly 200-300 during the day and ca 3000 during the night
    Opt_CR	=	optimset('TolFun',1); % [ppm] Numerical tolerance for internal CO2 computation
    if (LAI_roof > 0) 
        [rs_sun,rs_shd,Ci_sun,Ci_shd,~,~,~,~]=resistance_functions.Canopy_Resistence_An_Evolution(PAR_sun,PAR_shd,LAI_roof,...
            Kopt_roof,Knit_roof,Fsun,Fshd,Citm1_sun,Citm1_shd,...
            Catm_CO2,ra,rb,Troof_veg_tm1-273.15,Tatm-273.15,Pre/100,Ds_atm,...
            Psi_ltm1,Psi_sto_50_roof,Psi_sto_00_roof,...
            CT_roof,Vmax_roof,DSE_roof,Ha_roof,FI_roof,Oa,Do_roof,a1_roof,go_roof,e_rel_roof,...
            e_relN_roof,gmes_roof,rjv_roof,mSl_roof,Sl_roof,Opt_CR);
    else
        rs_sun=Inf;  rs_shd=Inf; Ci_sun=0; Ci_shd=0;
    end
end

Ci_sun_roof		=	Ci_sun;	% Ci = Leaf Interior  CO2 concentration [umolCO2/mol]
Ci_shd_roof		=	Ci_shd;

% Calcluate soil resistance
[~,~,~,Osat,Ohy,nVG,alpVG,Ks_Zs,L,Pe,O33,SPAR,~,~,...
	~,Rf_Zs,~,~,~,~,~,~,~,~,~,~]...
	=soil_functions.SoilParametersTotal(Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT,CASE_ROOT,0,ZR95,0,ZR50,0,ZRmax,Zs);

[r_soil,~,alp_soil]=resistance_functions.Soil_Resistence(Troof_veg-273.15,Pre/100,u_ref_und,ea,In_ground_tm1,Otm1(1),Ks_Zs(1),Osat(1),Ohy(1),L(1),Pe(1),O33(1),alpVG(1),nVG(1),SPAR);
% [r_soil,b_soil,alp_soil]=Soil_Resistence(Ts,Pre,Ws,ea,q_runon,O,Ks,Osat,Ohy,L,Pe,O33,alpVG,nVG,SPAR)

%% Calculation of turbulent fluxes
Hroof_imp			=	cp_atm*rho_atm*(Troof_imp-Tatm)/ra;	% sensible heat from impervious area [W/m^2]
Hroof_veg			=	cp_atm*rho_atm*(Troof_veg-Tatm)/(rb/(2*(LAI_roof+SAI_roof))+ra);	% sensible heat from vegetated area [W/m^2]

Eroof_imp_pot		=	rho_atm*(qsat_T_rimp-q_atm)/ra;	% Potential evaporation from ponding water on impervious area [kg/m^2 s]

Eroof_veg_pot		=	rho_atm*(qsat_T_rveg-q_atm)/(rb/((LAI_roof+SAI_roof)*dw_veg_roof)+ra);	% Potential evaporation from intercepted water on vegetation [kg/m^2 s]
%Eroof_ground_pot	=	rho_atm*(qsat_T_rveg-q_atm)/ra;	% Potential evaporation from ponding water on ground under vegetation [kg/m^2 s]
Eroof_soil_pot		=	rho_atm*alp_soil*(qsat_T_rveg-q_atm)/(ra+r_soil);	% Potential evaporation from first soil layer [kg/m^2 s]

TEroof_veg_sun_pot	=	rho_atm*(qsat_T_rveg-q_atm)/(rb/((LAI_roof)*(1-dw_veg_roof))+ra+rs_sun/((LAI_roof)*Fsun*(1-dw_veg_roof)));
TEroof_veg_shd_pot	=	rho_atm*(qsat_T_rveg-q_atm)/(rb/((LAI_roof)*(1-dw_veg_roof))+ra+rs_shd/((LAI_roof)*Fshd*(1-dw_veg_roof)));
TEroof_veg_pot		=	TEroof_veg_sun_pot+TEroof_veg_shd_pot;	% Total transpiration from vegetation  [kg/m^2 s]

% Condition that evapotranspiration does not exceed available water
Eroof_imp		=	min(Eroof_imp_pot,(In_imp_tm1/(1000*dth*3600)*row)); % Real max water evaporation from interception [kg/m^2 s]
Eroof_veg		=	min(Eroof_veg_pot,(In_veg_tm1/(1000*dth*3600)*row));
Eroof_ground	=	min(Eroof_soil_pot,(In_ground_tm1/(1000*dth*3600)*row));
Eroof_soil_pot	=	Eroof_soil_pot - Eroof_ground;

% Water available for Transpiration and Evaporation in a given time step
Vavail_tm1			=	(Vtm1./dth).*(row/3600/1000);	% Water volume in each soil layer [kg/m^2.s]
Eroof_soil			=	min(Eroof_soil_pot,Vavail_tm1(1));
Vavail_tm1(1)		=	Vavail_tm1(1)-Eroof_soil;
Vavail_plant_tm1	=	min(Vavail_tm1.*Rf_Zs,Exwat_tm1.*(row/3600/1000));	% Water volume for plant respiration in each soil layer [kg/m^2.s]
Vavail_plant_tm1 	=	sum(Vavail_plant_tm1,2);	% Cumulated water volume for plant transpiration [kg/m^2.s]
TEroof_veg			=	min(TEroof_veg_pot,Vavail_plant_tm1);

if fveg_roof==0
	Hroof_veg		=	0;
	Eroof_veg		=	0;
	Eroof_ground	=	0;
	Eroof_soil		=	0;
	TEroof_veg		=	0;
end

if fimp_roof==0
    Hroof_imp		=	0;
    Eroof_imp		=	0;
end

LEroof_imp			=	L_heat*Eroof_imp;		% Latent heat from ponding water on impervious area [W/m^2]
LEroof_veg			=	L_heat*Eroof_veg;		% Latent heat from intercepted water on vegetation [W/m^2]
LEroof_ground		=	L_heat*Eroof_ground;	% Latent heat from ponding water on ground under vegetation [W/m^2]
LEroof_soil			=	L_heat*Eroof_soil;		% Latent heat from first soil layer [W/m^2]
LTEroof_veg			=	L_heat*TEroof_veg;		% Latent heat from vegetation  [W/m^2]\

rap_L	=	NaN;
