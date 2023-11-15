function[rs_sun,rs_shd,Ci_sun,Ci_shd]=PrecalculateStomatalResistanceRoof(TempVec_ittm,MeteoData,HumidityAtm,...
    ParVegRoof,SoilPotW_ittm,CiCO2Leaf_ittm,PropOpticalRoof,ra,rb)

Troof_veg_tm1   =	TempVec_ittm.TRoofVeg;
Tatm			=	MeteoData.Tatm;
Pre				=	MeteoData.Pre;
ea				=	MeteoData.ea;
Catm_O2			=	MeteoData.Catm_O2;
Catm_CO2		=	MeteoData.Catm_CO2;
esat_Tatm		=	HumidityAtm.AtmVapourPreSat;
LAI_roof		=	ParVegRoof.LAI;
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

Psi_ltm1		=	SoilPotW_ittm.SoilPotWRoofVeg_L;
Ci_sun_tm1		=	CiCO2Leaf_ittm.CiCO2LeafRoofVegSun;
Ci_shd_tm1		=	CiCO2Leaf_ittm.CiCO2LeafRoofVegShd;


%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for stomata resistance
Citm1_sun	=	Ci_sun_tm1;		% Ci = Leaf Interior  CO2 concentration [umolCO2/mol]
Citm1_shd	=	Ci_shd_tm1;		% Ci = Leaf Interior  CO2 concentration [umolCO2/mol]
Ds_atm		=	esat_Tatm-ea;	% Ds = Vapor Pressure Deficit [Pa]
Oa			=	Catm_O2;		% Oa Intercellular Partial Pressure Oxygen [umolO2/mol]

SWRabs_dir		=	(1-PropOpticalRoof.aveg)*MeteoData.SW_dir;	% Absorbed direct shortwave radiation by vegetated roof [W/m^2]
SWRabs_diff		=	(1-PropOpticalRoof.aveg)*MeteoData.SW_diff;	% Absorbed diffuse shortwave radiation by vegetated roof [W/m^2]

% Partitioning of radiation into sunlit and shaded area
Fsun			=	(1.0 - exp(-Kopt_roof*(LAI_roof)))/(Kopt_roof*(LAI_roof));
Fsun(Fsun<0.01)	=	0; 
Fsun(Fsun>1)	=	1;
Fshd			=	1- Fsun;
PAR_sun			=	SWRabs_dir+Fsun*SWRabs_diff;	% [W/m^2] absorbed direct and diffuse shortwave radiation of the sunlit surface
PAR_shd			=	Fshd*SWRabs_diff;				% [W/m^2] absorbed direct and diffuse shortwave radiation of the shaded surface 


% Stomatal resistance,roughly 200-300 during the day and ca 3000 during the night
Opt_CR	=	optimset('TolFun',1); % [ppm] Numerical tolerance for internal CO2 computation
[rs_sun,rs_shd,Ci_sun,Ci_shd,~,~,~,~]=resistance_functions.Canopy_Resistence_An_Evolution(PAR_sun,PAR_shd,LAI_roof,...
    Kopt_roof,Knit_roof,Fsun,Fshd,Citm1_sun,Citm1_shd,...
    Catm_CO2,ra,rb,Troof_veg_tm1-273.15,Tatm-273.15,Pre/100,Ds_atm,...
    Psi_ltm1,Psi_sto_50_roof,Psi_sto_00_roof,...
    CT_roof,Vmax_roof,DSE_roof,Ha_roof,FI_roof,Oa,Do_roof,a1_roof,go_roof,e_rel_roof,...
    e_relN_roof,gmes_roof,rjv_roof,mSl_roof,Sl_roof,Opt_CR);






