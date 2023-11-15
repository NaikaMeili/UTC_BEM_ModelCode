function[q_runon_imp,In_imp,dIn_imp_dt,Lk_imp,...
	q_runon_veg,In_veg,dIn_veg_dt,...
	q_runon_ground,In_ground,dIn_ground_dt,dV_dt,f_ground,...
	V,O,OS,Lk,Psi_s,Exwat,Rd,TEroof_veg,Eroof_soil,Runoff,Runon_ittm,...
	WBalance_In_imp,WBalance_In_veg,WBalance_In_ground,WBalance_soil,...
	WBalance_imp_tot,WBalance_veg_tot,WBalance_tot]=...
	WaterRoof(Eroof_imp,Eroof_veg,Eroof_ground,Eroof_soil,TEroof_veg,...
	MeteoData,Int_ittm,Owater_ittm,Runon_ittm,FractionsRoof,ParSoilRoof,ParCalculation,ParVegRoof,Anthropogenic)


%% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rain				=	precipitaiton [mm]
% E_imp				=	evaporation from impervious surfaces [kg/m^2 s]
% E_veg				=	evaporation from vegetation [kg/m^2 s]
% E_ground			=	evaporation from ground under vegetation [kg/m^2 s]
% E_soil			=	evaporation from soil under vegetation [kg/m^2 s]
% TE_veg			=	Transpiration from vegetation [kg/m^2 s]
% In_imp_tm1		=	Interception from previous time step [mm]
% In_veg_tm1		=	Interception from previous time step [mm]
% In_ground_tm1		=	Interception from previous time step [mm]
% Otm1				=	Soil moisture / water content in the different soil layers [-]
% OStm1				=	[] Water Content for evaporation of previous time step
% In_max_imp		=	Maximum interception capacity of urban area [mm]
% In_max_ground		=	Maximum interception capacity of ground [mm]
% Sp_In				=	specific water retained by a vegetated surface [mm m^2 VEG area m^-2 plant area]
% K_imp				=	Hydraulic conductivity of impervious area [mm/h]
% dth				=	time step of computation [h]
% row				=	Water desity [kg/m^3]
% LAI				=	Leaf area index [-]
% SAI				=	Stem area index [-]
% Zs				=	soil layer discretization [mm]
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

%% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% q_runon_imp		=	Runoff [mm/h]
% In_imp			=	Interception [mm]
% dIn_imp_dt		=	Change in interception [mm]
% Lk_imp			=	Leakage from impervious area [mm/h]
% q_runon_veg		=	Runoff [mm/h]
% In_veg			=	Interception [mm]
% dIn_veg_dt		=	Change in interception [mm]
% q_runon_ground	=	Runoff [mm/h]
% In_ground			=	Interception [mm]
% dIn_ground_dt		=	Change in interception [mm]
% f_ground			=	Infiltration into first soil layer [mm/h]
% V					=	Volume of water in each soil layer [mm]. Never includes the residual water content Ohy
% O					=	Soil moisture in each soil layer [-]. Always includes the residual water content Ohy
% OS				=	[] Water Content for evaporation
% Lk				=	Leakage at bedrock [mm/h]
% Rd				=	saturation excess runoff / Dunne Runoff [mm] 
% Psi_s				=	Soil Water Potential for Vegetation [MPa] 
% Exwat				=	Max extractable water for Vegetation [mm m2 / m2 ground h ] 
% TE_veg			=	Transpiraton [kg/m^2 s]
% E_soil			=	Evaporation from first soil layer [kg/m^2 s]
% WBalance_In_imp	=	Water balance interception on impervious area [mm]
% WBalance_In_veg	=	Water balance of interception on vegetation [mm]
% WBalance_In_ground	= Water balance of ponding on the ground [mm]
% WBalance_soil		=	Water balance in the soil  [mm]

Rain			=	MeteoData.Rain;
In_imp_tm1		=	Int_ittm.IntRoofImp;
In_veg_tm1		=	Int_ittm.IntRoofVegPlant;
In_ground_tm1	=	Int_ittm.IntRoofVegGround;
Otm1			=	Owater_ittm.OwRoofSoilVeg;
Runon_tm1		=	Runon_ittm.RunonRoofTot;
Per_runoff		=	FractionsRoof.Per_runoff;
fveg			=	FractionsRoof.fveg;
fimp			=	FractionsRoof.fimp;
In_max_imp		=	ParSoilRoof.In_max_imp;
In_max_ground	=	ParSoilRoof.In_max_ground;
K_imp			=	ParSoilRoof.Kimp;
dth				=	ParCalculation.dth;
row				=	ParCalculation.row;
Sp_In			=	ParSoilRoof.Sp_In;
Zs				=	ParSoilRoof.Zs;
Pcla			=	ParSoilRoof.Pcla;
Psan			=	ParSoilRoof.Psan;
Porg			=	ParSoilRoof.Porg;
Kfc				=	ParSoilRoof.Kfc;
Phy				=	ParSoilRoof.Phy;
SPAR			=	ParSoilRoof.SPAR;
Kbot			=	ParSoilRoof.Kbot;
LAI				=	ParVegRoof.LAI;
SAI				=	ParVegRoof.SAI;
Rrootl			=	ParVegRoof.Rrootl;
PsiL50			=	ParVegRoof.PsiL50;
PsiX50			=	ParVegRoof.PsiX50;
CASE_ROOT		=	ParVegRoof.CASE_ROOT;
ZR95			=	ParVegRoof.ZR95;
ZR50			=	ParVegRoof.ZR50;
ZRmax			=	ParVegRoof.ZRmax;


%% Water flux compuations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[q_runon_imp,In_imp,dIn_imp_dt,Lk_imp,WBalance_In_imp]=...
	water_functions.Water_Impervious(Rain,Runon_tm1,Eroof_imp,In_imp_tm1,dth,row,In_max_imp,K_imp);

[q_runon_veg,In_veg,dIn_veg_dt,WBalance_In_veg]=...
	water_functions.Water_Vegetation(Rain,Eroof_veg,In_veg_tm1,Sp_In,LAI,SAI,row,dth);

[q_runon_ground,In_ground,dIn_ground_dt,f_ground,WBalance_In_ground]=...
	water_functions.Water_Ground(q_runon_veg+Anthropogenic.Waterf_roof,Runon_tm1,...
	Eroof_ground,Otm1,In_ground_tm1,In_max_ground,...
	Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,CASE_ROOT,CASE_ROOT,0,ZR95,0,ZR50,0,ZRmax,Zs,dth,row);

[V,O,OS,Lk,~,Psi_s,~,Exwat,Rd,TEroof_veg,~,Eroof_soil,dV_dt,WBalance_soil,Psi_soil,Ko]...
	=water_functions.Water_Soil(Otm1,f_ground,0,TEroof_veg,Eroof_soil,0.*Zs,...
	dth,Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,CASE_ROOT,CASE_ROOT,0,ZR95,0,ZR50,0,ZRmax,...
	0,Rrootl,0,PsiL50,0,PsiX50,Zs,row);

Runoff				=	Per_runoff*(fimp*q_runon_imp + fveg*(q_runon_ground + Rd));			% [mm/dth]
Runon_ittm				=	(1 - Per_runoff)*(fimp*q_runon_imp + fveg*(q_runon_ground + Rd));	% [mm/dth]

% Calculate water balance for check
WBalance_imp_tot	=	Rain + Runon_tm1 - Eroof_imp*dth*3600*1000/row - q_runon_imp - Lk_imp*dth - dIn_imp_dt;	% [mm/dth]

WBalance_veg_tot	=	Rain + Runon_tm1 + Anthropogenic.Waterf_roof - (Eroof_veg + Eroof_ground + Eroof_soil + TEroof_veg)*dth*3600*1000/row - Lk*dth - q_runon_ground - Rd...
						- dIn_veg_dt - dIn_ground_dt - dV_dt;	% [mm/dth]
					
E_tot				=	(fimp*Eroof_imp + fveg*(Eroof_veg + Eroof_ground + Eroof_soil + TEroof_veg))*3600*1000/row;	% [mm/h]
Leak_tot			=	fimp*Lk_imp + fveg*Lk;	% [mm/h]
Storage_tot			=	fimp*dIn_imp_dt + fveg*(dIn_veg_dt + dIn_ground_dt + dV_dt);	% [mm/dth]
					
WBalance_tot		=	Rain + Runon_tm1 + fveg*Anthropogenic.Waterf_roof - E_tot*dth - Leak_tot*dth - Runoff - Runon_ittm - Storage_tot;	% [mm/dth]

end
