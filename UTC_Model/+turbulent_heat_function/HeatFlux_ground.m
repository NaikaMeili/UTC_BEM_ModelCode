function[Himp,Hbare,Hveg,Htree,...
	Eimp,Ebare_pond,Ebare_soil,Eveg_int,...
	Eveg_pond,Eveg_soil,TEveg,Etree_int,TEtree,...
	Ebare,Eveg,Etree,...
	LEimp,LEbare_pond,LEbare_soil,LEveg_int,...
	LEveg_pond,LEveg_soil,LTEveg,LEtree_int,LTEtree,...
	LEbare,LEveg,LEtree,...
	Ci_sun_H,Ci_shd_H,Ci_sun_L,Ci_shd_L,...
	rap_can,rap_Htree_In,rb_H,rb_L,...
	r_soil_bare,r_soil_veg,alp_soil_bare,alp_soil_veg,...
	rs_sun_L,rs_shd_L,rs_sun_H,rs_shd_H,...
	u_Hcan,u_Zref_und,Fsun_L,Fshd_L,dw_L]...
    =HeatFlux_ground(TemperatureC,TempVec_ittm,MeteoData,Gemeotry_m,geometry,FractionsGround,ParTree,...
	ParVegGround,ParVegTree,ParSoilGround,SoilPotW_ittm,Owater_ittm,Vwater_ittm,ExWater_ittm,Int_ittm,CiCO2Leaf_ittm,...
	ParInterceptionTree,ParCalculation,SWRdir_abs_tree,SWRdiff_abs_tree,SWRdir_abs_groundveg,SWRdiff_abs_groundveg,...
    RESPreCalc,rsGroundPreCalc,rsTreePreCalc)


%% Parameter specification
Timp			=	TemperatureC(1,1);
Tbare			=	TemperatureC(1,2);
Tveg			=	TemperatureC(1,3);
Ttree			=	TemperatureC(1,6);
Tcanyon			=	TemperatureC(1,9);
Twsun			=	TemperatureC(1,4);
Twshade			=	TemperatureC(1,5);
qcanyon			=	TemperatureC(1,10);
Tveg_tm1		=	TempVec_ittm.TGroundVeg;
Ttree_tm1		=	TempVec_ittm.TTree;
Tcanyon_tm1     =   TempVec_ittm.TCanyon;
Tatm			=	MeteoData.Tatm;
Pre				=	MeteoData.Pre;
ea				=	MeteoData.ea;
Zatm			=	MeteoData.Zatm;
Uatm			=	MeteoData.Uatm;
Catm_O2			=	MeteoData.Catm_O2;
Catm_CO2		=	MeteoData.Catm_CO2;
H				=	Gemeotry_m.Height_canyon;
W				=	Gemeotry_m.Width_canyon;
w				=	geometry.wcanyon;
Wroof			=	Gemeotry_m.Width_roof;
Htree			=	Gemeotry_m.Height_tree;
R_tree			=	Gemeotry_m.Radius_tree;
Hcan_max		=	Gemeotry_m.Hcan_max;
Hcan_std		=	Gemeotry_m.Hcan_std;
wroof_norm		=	geometry.wroof_norm;	
rad_tree		=	geometry.radius_tree;
fgveg			=	FractionsGround.fveg;
fgbare			=	FractionsGround.fbare;
fgimp			=	FractionsGround.fimp;
trees			=	ParTree.trees;
LAI_L			=	ParVegGround.LAI;
SAI_L			=	ParVegGround.SAI;
LAI_H			=	ParVegTree.LAI;
SAI_H			=	ParVegTree.SAI;
hc_L			=	ParVegGround.hc;
hc_H			=	Gemeotry_m.Height_tree;
d_leaf_L		=	ParVegGround.d_leaf;
d_leaf_H		=	ParVegTree.d_leaf;
Kopt_H			=	ParVegTree.Kopt;
Kopt_L			=	ParVegGround.Kopt;
Knit_H			=	ParVegTree.Knit;
Knit_L			=	ParVegGround.Knit;
Psi_sto_50_H	=	ParVegTree.Psi_sto_50;
Psi_sto_50_L	=	ParVegGround.Psi_sto_50;
Psi_sto_00_H	=	ParVegTree.Psi_sto_00;
Psi_sto_00_L	=	ParVegGround.Psi_sto_00;
CT_H			=	ParVegTree.CT;
CT_L			=	ParVegGround.CT;
Vmax_H			=	ParVegTree.Vmax;
Vmax_L			=	ParVegGround.Vmax;
DSE_H			=	ParVegTree.DSE;
DSE_L			=	ParVegGround.DSE;
Ha_H			=	ParVegTree.Ha;
Ha_L			=	ParVegGround.Ha;
FI_H			=	ParVegTree.FI;
FI_L			=	ParVegGround.FI;
Do_H			=	ParVegTree.Do;
Do_L			=	ParVegGround.Do;
a1_H			=	ParVegTree.a1;
a1_L			=	ParVegGround.a1;
go_H			=	ParVegTree.go;
go_L			=	ParVegGround.go;
e_rel_H			=	ParVegTree.e_rel;
e_rel_L			=	ParVegGround.e_rel;
e_relN_H		=	ParVegTree.e_relN;
e_relN_L		=	ParVegGround.e_relN;
gmes_H			=	ParVegTree.gmes;
gmes_L			=	ParVegGround.gmes;
rjv_H			=	ParVegTree.rjv;
rjv_L			=	ParVegGround.rjv;
mSl_H			=	ParVegTree.mSl;
mSl_L			=	ParVegGround.mSl;
Sl_H			=	ParVegTree.Sl;
Sl_L			=	ParVegGround.Sl;
SPARTREE		=	ParVegTree.SPARTREE;
Pcla			=	ParSoilGround.Pcla;
Psan			=	ParSoilGround.Psan;
Porg			=	ParSoilGround.Porg;
Kfc				=	ParSoilGround.Kfc;
Phy				=	ParSoilGround.Phy;
SPAR			=	ParSoilGround.SPAR;
Kbot			=	ParSoilGround.Kbot;
CASE_ROOT_H		=	ParVegTree.CASE_ROOT;
CASE_ROOT_L		=	ParVegGround.CASE_ROOT;
ZR95_H			=	ParVegTree.ZR95;
ZR95_L			=	ParVegGround.ZR95;
ZR50_H			=	ParVegTree.ZR50;
ZR50_L			=	ParVegGround.ZR50;
ZRmax_H			=	ParVegTree.ZRmax;
ZRmax_L			=	ParVegGround.ZRmax;
Zs				=	ParSoilGround.Zs;
Psi_L_tm1		=	SoilPotW_ittm.SoilPotWGroundVeg_L;
Psi_H_tm1		=	SoilPotW_ittm.SoilPotWGroundTot_H;
Otm1Imp			=	Owater_ittm.OwGroundSoilImp;
Otm1Bare		=	Owater_ittm.OwGroundSoilBare;
Otm1Veg			=	Owater_ittm.OwGroundSoilVeg;
Vtm1Imp			=	Vwater_ittm.VGroundSoilImp;
Vtm1Bare		=	Vwater_ittm.VGroundSoilBare;
Vtm1Veg			=	Vwater_ittm.VGroundSoilVeg;
dth				=	ParCalculation.dth;
row				=	ParCalculation.row;
ExwatImp_tm1_H	=	ExWater_ittm.ExWaterGroundImp_H;
ExwatBare_tm1_H	=	ExWater_ittm.ExWaterGroundBare_H;
ExwatVeg_tm1_H	=	ExWater_ittm.ExWaterGroundVeg_H;
ExwatVeg_tm1_L	=	ExWater_ittm.ExWaterGroundVeg_L;
Sp_In_H			=	ParInterceptionTree.Sp_In;
Sp_In_L			=	ParSoilGround.Sp_In;
In_imp_tm1		=	Int_ittm.IntGroundImp;
In_bare_tm1		=	Int_ittm.IntGroundBare;
In_veg_tm1		=	Int_ittm.IntGroundVegPlant;
In_underveg_tm1	=	Int_ittm.IntGroundVegGround;
In_tree_tm1		=	Int_ittm.IntTree;
Ci_sun_H_tm1	=	CiCO2Leaf_ittm.CiCO2LeafTreeSun;
Ci_shd_H_tm1	=	CiCO2Leaf_ittm.CiCO2LeafTreeShd;
Ci_sun_L_tm1	=	CiCO2Leaf_ittm.CiCO2LeafGroundVegSun;
Ci_shd_L_tm1	=	CiCO2Leaf_ittm.CiCO2LeafGroundVegShd;


%% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hground_imp			=	sensible heat of impervious ground fraction (W/m^2)
% Hground_bare			=	sensible heat of bare ground fraction (W/m^2)
% Hground_veg			=	sensible heat of vegetated ground fraction (W/m^2)
% Eground_imp			=	latent heat of impervious ground fraction (x)
% Eground_bare			=	latent heat of bare ground fraction (x)
% Eground_veg_evapo		=	latent heat of evaporation of vegetated ground fraction (W/m^2)
% Eground_veg_transp	=	latent heat of transpiration of vegetated ground fraction (W/m^2)

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H						=	Canyon (building) height (m)
% w						=	ground width (-)
% radius_a				=	normalized tree radius (-)
% fgveg					=	Partitioning ground vegetation(-)
% fgbare				=	Partitioning ground bare(-
% fgimp					=	Partitioning ground impervious(-)
% LAI_tree				=	Leaf area index for the trees (-)
% LAI_ground			=	Leaf area index for the ground vegetation (-)
% SAI_ground			=	Stem area index for the ground vegetation (-)
% SAI_tree				=	Stem area index for the trees (-)
% hc_tree				=	canopy height trees [m]
% hc_ground				=	canopy height ground vegetation [m]
% d_leaf_tree			=	Leaf dimension of trees [cm]
% d_leaf_ground			=	Leaf dimension of ground vegetation [cm]
% Pre					=	air pressure [Pa]. Carefull, used to be [hPa - mbar] in T&C
% U_top_canyon			=	horizontal wind speed at the top of the canyon [m/s]
% cp_atm				=	specific heat air  [J/kg K]
% rho_atm				=	dry air density at atmosphere [kg/m^3]
% Tground_imp			=	Temperature of the ground impervious [K]
% Tground_bare			=	Temperature of the ground bare [K]
% Tground_veg			=	Temperature of the ground vegetated [K]
% Ttree					=	Temperature of the tree [K]
% T_canyon				=	Temperature of canyon air [K]
% Tground_tree			=	average projected ground temperature [K]
% qsat_T_gimp			=	Saturated specific humidity at Tground_imp []
% qsat_T_gbare			=	Saturated specific humidity at Tground_bare []
% qsat_T_gveg			=	Saturated specific humidity at Tground_veg []
% q_canyon				=	Specific humidity of canyon air []

% Structural parameters
Cimp			=	(fgimp>0);
Cbare			=	(fgbare>0);
Cveg			=	(fgveg>0);
Ctree			=	(trees==1);
Ctree025		=	(4*rad_tree>= 0.25*w);

if isnan(rad_tree)==1
	rad_tree	=	0;
end

% hc_ground		=	Cveg*hc_ground;
% LAI_ground		=	Cveg*LAI_ground;
% d_leaf_ground	=	Cveg*d_leaf_ground;
% Ht				=	Ctree*Ht;
% LAI_tree		=	Ctree*LAI_tree; 
% d_leaf_tree		=	Ctree*d_leaf_tree; 
% radius_a		=	Ctree*radius_a;

%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cp_atm			=	1005+(((Tatm-273.15)+23.15)^2)/3364;		% specific heat air  [J/kg K]
rho_atm			=	(Pre/(287.04*Tatm))*(1-(ea/Pre)*(1-0.622));	% dry air density at atmosphere [kg/m^3]
q_atm			=	0.622*ea/(Pre-0.378*ea);					% Specifc humidity of air at reference height []
L_heat			=	1000*(2501.3 - 2.361*(Tatm-273.15));		% Latent heat vaporization/condensation [J/kg]

% Average temperatures
Tgroundtree		=	(fgveg*Tveg+fgbare*Tbare+fgimp*Timp+Ctree*Ttree*(4*rad_tree))/(fgveg+fgbare+fgimp+(Ctree*4*rad_tree));    % Average temperature of the ground [K]

Tsurf		=	(fgveg*Tveg+fgbare*Tbare+fgimp*Timp+Ctree*Ttree*(4*rad_tree)+H/W*Twsun+H/W*Twshade)/(fgveg+fgbare+fgimp+(Ctree*4*rad_tree)+2*H/W);    % Average temperature of the ground [K]

Tground			=	fgveg*Tveg+fgbare*Tbare+fgimp*Timp;    % Average temperature of the ground [K]
Tgroundsoil		=	(fgveg*Tveg+fgbare*Tbare)/(fgveg+fgbare);    % Average temperature of the ground [K]

% Vapor pressure saturation and specific humidity at esat
esat_T_imp		=	611*exp(17.27*(Timp-273.16)/(237.3+(Timp-273.16)));	% vapor pressure saturation at Tground_imp [Pa]
qsat_T_imp		=	(0.622*esat_T_imp)/(Pre-0.378*esat_T_imp);						% Saturated specific humidity at Tground_imp []
esat_T_bare		=	611*exp(17.27*(Tbare-273.16)/(237.3+(Tbare-273.16)));	% vapor pressure saturation at Tground_bare [Pa]
qsat_T_bare		=	(0.622*esat_T_bare)/(Pre-0.378*esat_T_bare);						% Saturated specific humidity at Tground_bare []
esat_T_veg		=	611*exp(17.27*(Tveg-273.16)/(237.3+(Tveg-273.16)));	% vapor pressure saturation at Tground_veg [Pa]
qsat_T_veg		=	(0.622*esat_T_veg)/(Pre-0.378*esat_T_veg);						% Saturated specific humidity at Tground_veg []
esat_T_tree		=	611*exp(17.27*(Ttree-273.16)/(237.3+(Ttree-273.16)));				% vapor pressure saturation at Ttree [Pa]
qsat_T_tree		=	(0.622*esat_T_tree)/(Pre-0.378*esat_T_tree);						% Saturated specific humidity at Ttree []
esat_T_canyon	=	611*exp(17.27*(Tcanyon-273.16)/(237.3+(Tcanyon-273.16)));			% vapor pressure saturation at T_canyon [Pa]
qsat_T_canyon	=	(0.622*esat_T_canyon)/(Pre-0.378*esat_T_canyon);					% Saturated specific humidity at T_canyon []
e_T_canyon		=	qcanyon*Pre/(0.622+0.378*qcanyon);
rel_hum_canyon	=	e_T_canyon/esat_T_canyon;
esat_Timpbarevegtree	=	611*exp(17.27*(Tgroundtree-273.16)/(237.3+(Tgroundtree-273.16)));				% vapor pressure saturation at Ttree [Pa]

% Parameters for stomata resistance
Citm1_sun_H		=	Ctree.*Ci_sun_H_tm1;% Ci = Leaf Interior  CO2 concentration [umolCO2/mol]
Citm1_shd_H		=	Ctree.*Ci_shd_H_tm1;% Ci = Leaf Interior  CO2 concentration [umolCO2/mol]
Citm1_sun_L		=	Cveg.*Ci_sun_L_tm1;	% Ci = Leaf Interior  CO2 concentration [umolCO2/mol]
Citm1_shd_L		=	Cveg.*Ci_shd_L_tm1;	% Ci = Leaf Interior  CO2 concentration [umolCO2/mol]

Ds_canyon		=	esat_T_canyon - e_T_canyon;		% Ds = Vapor Pressure Deficit [Pa]
Oa				=	Catm_O2;						% Oa Intercellular Partial Pressure Oxygen [umolO2/mol]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensible heat flux ground
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partitioning of radiation into sunlit and shaded area
Fsun_H				=	Ctree*((1.0 - exp(-Kopt_H*(LAI_H)))/(Kopt_H*(LAI_H)));
Fsun_H(Fsun_H<0.01)	=	0; 
Fsun_H(Fsun_H>1)	=	1;
Fshd_H				=	Ctree*(1- Fsun_H);
PAR_sun_H			=	Ctree*(SWRdir_abs_tree + Fsun_H*SWRdiff_abs_tree);		% [W/m^2] absorbed direct and diffuse shortwave radiation of the sunlit surface
PAR_shd_H			=	Ctree*(Fshd_H*SWRdiff_abs_tree);						% [W/m^2] absorbed direct and diffuse shortwave radiation of the shaded surface 

Fsun_L				=	Cveg*((1.0 - exp(-Kopt_L*(LAI_L)))/(Kopt_L*(LAI_L)));
Fsun_L(Fsun_L<0.01)	=	0; 
Fsun_L(Fsun_L>1)	=	1;
Fshd_L				=	Cveg*(1- Fsun_L);
PAR_sun_L			=	Cveg*(SWRdir_abs_groundveg+Fsun_L*SWRdiff_abs_groundveg);	% [W/m^2] absorbed direct and diffuse shortwave radiation of the sunlit surface
PAR_shd_L			=	Cveg*(Fshd_L*SWRdiff_abs_groundveg);						% [W/m^2] absorbed direct and diffuse shortwave radiation of the shaded surface 

% Fraction of vegetation covered by intercepted water (Deardorff(1978))
In_max_H			=	Ctree*(Sp_In_H*(LAI_H+SAI_H));
dw_H				=	Ctree*(min(1,(In_tree_tm1/In_max_H)^(2/3)));				% Fraction of vegetation covered by intercepted water
In_max_L			=	Cveg*(Sp_In_L*(LAI_L+SAI_L));
dw_L				=	Cveg*(min(1,(In_veg_tm1/In_max_L)^(2/3)));	% Fraction of vegetation covered by intercepted water

% Calcualtion of structural parameters and wind profile in the city
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[zom,zoh,zom_ground,zoh_ground,disp_h,zom_H,zom_L,zoh_H,zoh_L,d_H,d_L,zom_other]...
	=resistance_functions.Urban_roughness(Ctree*hc_H,Cveg*hc_L,Cbare,Cimp,0);
%[zom,zoh,disp_h,zom_H,zom_L,zoh_H,zoh_L,d_H,d_L,zom_other]=Urban_roughness(hc_H,hc_L,Csoil,Croad,Croof)

if Ctree>0
[dcan,zomcan,u_Hcan,u_tree,w_tree,alpha,RoughnessParameter]=resistance_functions.WindProfile_Canyon(...
	H,Htree,R_tree,W,Wroof,Kopt_H,LAI_H,Zatm,Uatm,Ctree*hc_H,trees,1.5,zom_ground,Hcan_max,Hcan_std);
else
	u_tree=0.1;
end

if Cveg>0
[dcan,zomcan,u_Hcan,u_Lveg,w_Lveg,alpha,RoughnessParameter]=resistance_functions.WindProfile_Canyon(...
	H,Htree,R_tree,W,Wroof,Kopt_H,LAI_H,Zatm,Uatm,Cveg*hc_L,trees,1.5,zom_ground,Hcan_max,Hcan_std);
else
	u_Lveg=0.1;
end

[dcan,zomcan,u_Hcan,u_Zref_und,w_Zref_und,alpha,RoughnessParameter]=resistance_functions.WindProfile_Canyon(...
	H,Htree,R_tree,W,Wroof,Kopt_H,LAI_H,Zatm,Uatm,1.5,trees,1.5,zom_ground,Hcan_max,Hcan_std);
% [dcan,zomcan,u_Hcan,u_Zp,w_Zp]=...
% 	WindProfile_Canyon(Hcan,Htree,R_tree,Wcan,Wroof,Kopt,LAI_t,Zatm,uatm,Zp,trees,Zref_und,zom_und)

%% Calculation of aerodynamic and stomata resistances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rap_can,rap_Htree,rap_Htree_In,rap_2m,rap_2m_In,rap_Zp3,rap_Zp3_In,u_Hcan,u_Htree,u_Zp2,u_Zp3,uref_und,alpha]...
	=resistance_functions.InCanyonAerodynamicResistance(Uatm,Zatm,Tcanyon-273.15,Tsurf-273.15,...
	Hcan_max,H,dcan,zomcan,1.5,zom_ground,hc_H,2,2,Pre,e_T_canyon,RoughnessParameter);
% [rap_can,rap_Zp1,rap_Zp1_In,rap_Zp2,rap_Zp2_In,rap_Zp3,rap_Zp3_In,u_Hcan,u_Zp1,u_Zp2,u_Zp3,uref_und]...
% 	=InCanyonAerodynamicResistance(uatm,Zatm,Ta,Ts,hcan,dcan,zomcan,Zref_und,zom_und,Zp1,Zp2,Zp3)

Opt_CR	=	optimset('TolFun',1); % [ppm] Numerical tolerance for internal CO2 computation

if Ctree==1 && LAI_H>0
	[rb_H]=resistance_functions.Leaf_BR(u_tree,Ttree-273.15,Tcanyon-273.15,d_leaf_H,alpha);
	% [rb]=Leaf_BR(u_hc,Ts,Ta,d_leaf,alpha)

    if RESPreCalc==1
        rs_sun_H = rsTreePreCalc.rs_sun_H;  
        rs_shd_H = rsTreePreCalc.rs_shd_H; 
        Ci_sun_H = rsTreePreCalc.Ci_sun_H; 
        Ci_shd_H = rsTreePreCalc.Ci_shd_H;
    else
	    [rs_sun_H,rs_shd_H,Ci_sun_H,Ci_shd_H,~,~,~,~]=resistance_functions.Canopy_Resistence_An_Evolution(PAR_sun_H,PAR_shd_H,LAI_H,...
		    Kopt_H,Knit_H,Fsun_H,Fshd_H,Citm1_sun_H,Citm1_shd_H,...
		    Catm_CO2,rap_Htree_In,rb_H,Ttree_tm1-273.15,Tcanyon_tm1-273.15,Pre/100,Ds_canyon,...
		    Psi_H_tm1,Psi_sto_50_H,Psi_sto_00_H,...
		    CT_H,Vmax_H,DSE_H,Ha_H,FI_H,Oa,Do_H,a1_H,go_H,e_rel_H,...
		    e_relN_H,gmes_H,rjv_H,mSl_H,Sl_H,Opt_CR);
    end
else
	rb_H=Inf;	rs_sun_H=Inf;  rs_shd_H=Inf; Ci_sun_H=0; Ci_shd_H=0;
end

if Cveg==1 && LAI_L>0
	[rb_L]=resistance_functions.Leaf_BR(u_Lveg,Tveg-273.15,Tcanyon-273.15,d_leaf_L,alpha);
	% [rb]=Leaf_BR(u_hc,Ts,Ta,d_leaf,alpha)
	
% 	[~,rap_L,~,~,~]=resistance_functions.Undercanopy_Leaf_Resistence2(u_Zref_und,Tcanyon-273.15,Tground-273.15,1,0,Cveg*hc_L,0,LAI_L,0,d_leaf_L,...
%     1.5,d_L,zom_L,zoh_ground,0,0,0,d_L,zom_L);

    if RESPreCalc==1
        rs_sun_L = rsGroundPreCalc.rs_sun_L;  
        rs_shd_L = rsGroundPreCalc.rs_shd_L; 
        Ci_sun_L = rsGroundPreCalc.Ci_sun_L; 
        Ci_shd_L = rsGroundPreCalc.Ci_shd_L;
    else
	    [rs_sun_L,rs_shd_L,Ci_sun_L,Ci_shd_L,~,~,~,~]=resistance_functions.Canopy_Resistence_An_Evolution(PAR_sun_L,PAR_shd_L,LAI_L,...
		    Kopt_L,Knit_L,Fsun_L,Fshd_L,Citm1_sun_L,Citm1_shd_L,...
		    Catm_CO2,rap_can,rb_L,Tveg_tm1-273.15,Tcanyon_tm1-273.15,Pre/100,Ds_canyon,...
		    Psi_L_tm1,Psi_sto_50_L,Psi_sto_00_L,...
		    CT_L,Vmax_L,DSE_L,Ha_L,FI_L,Oa,Do_L,a1_L,go_L,e_rel_L,...
		    e_relN_L,gmes_L,rjv_L,mSl_L,Sl_L,Opt_CR);
    end
else
	rb_L=Inf;	rs_sun_L=Inf;  rs_shd_L=Inf; Ci_sun_L=0; Ci_shd_L=0;
end

%% Calculation of soil resistance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,~,Osat,Ohy,nVG,alpVG,Ks_Zs,L,Pe,O33,SPAR,~,~,RfH_Zs,RfL_Zs,~,~,~,~,~,~,~,~,~,~]...
	=soil_functions.SoilParametersTotal(Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,Zs);

[r_soil_bare,~,alp_soil_bare]=resistance_functions.Soil_Resistence...
	(Tbare-273.15,Pre/100,u_Zref_und,e_T_canyon,In_bare_tm1,Otm1Bare(1),...
	Ks_Zs(1),Osat(1),Ohy(1),L(1),Pe(1),O33(1),alpVG(1),nVG(1),SPAR);

[r_soil_veg,~,alp_soil_veg]=resistance_functions.Soil_Resistence...
	(Tveg-273.15,Pre/100,u_Zref_und,e_T_canyon,In_bare_tm1,Otm1Veg(1),...
	Ks_Zs(1),Osat(1),Ohy(1),L(1),Pe(1),O33(1),alpVG(1),nVG(1),SPAR);
% [r_soil,b_soil,alp_soil]=Soil_Resistence(Ts,Pre,Ws,ea,q_runon,O,Ks,Osat,Ohy,L,Pe,O33,alpVG,nVG,SPAR)


%% Calculation of turbulent fluxes
% Turbulent fluxes 
Himp			=	Cimp*(cp_atm*rho_atm*(Timp-Tcanyon)/rap_can);
Hbare			=	Cbare*(cp_atm*rho_atm*(Tbare-Tcanyon)/rap_can);
Hveg			=	Cveg*(cp_atm*rho_atm*(Tveg-Tcanyon)/(rb_L/(2*(LAI_L+SAI_L))+rap_can));
Htree			=	Ctree*(cp_atm*rho_atm*(Ttree-Tcanyon)/(rb_H/(2*(LAI_H+SAI_H))+rap_Htree_In));

Eimp_pot		=	Cimp*(rho_atm*(qsat_T_imp-qcanyon)./rap_can);

%Ebare_pond_pot	=	Cbare*(rho_atm*(qsat_T_bare-qcanyon)./rap_can);
Ebare_soil_pot	=	Cbare*(rho_atm*(alp_soil_bare*qsat_T_bare-qcanyon)./(rap_can+r_soil_bare));

Eveg_int_pot	=	Cveg*(rho_atm*(qsat_T_veg-qcanyon)./(rb_L/((LAI_L+SAI_L)*dw_L)+rap_can));
% Eveg_pond_pot	=	Cveg*(rho_atm*(qsat_T_veg-qcanyon)./rap_can);
Eveg_soil_pot	=	Cveg*(rho_atm*(alp_soil_veg*qsat_T_veg-qcanyon)./(rap_can+r_soil_veg));

TEveg_sun_pot	=	Cveg*(rho_atm*(qsat_T_veg-qcanyon)./(rb_L/((LAI_L)*Fsun_L*(1-dw_L))+rap_can+rs_sun_L/((LAI_L)*Fsun_L*(1-dw_L))));
TEveg_shd_pot	=	Cveg*(rho_atm*(qsat_T_veg-qcanyon)./(rb_L/((LAI_L)*Fshd_L*(1-dw_L))+rap_can+rs_shd_L/((LAI_L)*Fshd_L*(1-dw_L))));
TEveg_pot		=	Cveg*(TEveg_sun_pot+TEveg_shd_pot);

Etree_int_pot	=	Ctree*(rho_atm*(qsat_T_tree-qcanyon)./(rb_H/((LAI_H+SAI_H)*dw_H)+rap_Htree_In));

TEtree_sun_pot	=	Ctree*(rho_atm*(qsat_T_tree-qcanyon)./(rb_H/((LAI_H)*Fsun_H*(1-dw_H))+rap_Htree_In+rs_sun_H/((LAI_H)*Fsun_H*(1-dw_H))));
TEtree_shd_pot	=	Ctree*(rho_atm*(qsat_T_tree-qcanyon)./(rb_H/((LAI_H)*Fshd_H*(1-dw_H))+rap_Htree_In+rs_shd_H/((LAI_H)*Fsun_H*(1-dw_H))));
TEtree_pot		=	Ctree*(TEtree_sun_pot+TEtree_shd_pot);


% Condition that evapotranspiration does not exceed available water
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Water limitations of interception and ponding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eimp			=	min(Eimp_pot,(In_imp_tm1/(1000*dth*3600)*row));				% Real max water evaporation from interception [kg/m^2 s]
Ebare_pond		=	min(Ebare_soil_pot,(In_bare_tm1/(1000*dth*3600)*row));		% Real max water evaporation from interception [kg/m^2 s]
Ebare_soil_pot	=	Ebare_soil_pot-Ebare_pond;
Eveg_int		=	min(Eveg_int_pot,(In_veg_tm1/(1000*dth*3600)*row));			% Real max water evaporation from interception [kg/m^2 s]
Eveg_pond		=	min(Eveg_soil_pot,(In_underveg_tm1/(1000*dth*3600)*row));	% Real max water evaporation from interception [kg/m^2 s]
Eveg_soil_pot	=	Eveg_soil_pot-Eveg_pond;
Etree_int		=	min(Etree_int_pot,(In_tree_tm1/(1000*dth*3600)*row));		% Real max water evaporation from interception [kg/m^2 s]

% Water limitation to soil evaporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VavailImp_tm1		=	(Vtm1Imp./dth).*(row/3600/1000);	% Water volume in each soil layer [kg/m^2.s]
VavailBare_tm1		=	(Vtm1Bare./dth).*(row/3600/1000);	% Water volume in each soil layer [kg/m^2.s]
VavailVeg_tm1		=	(Vtm1Veg./dth).*(row/3600/1000);	% Water volume in each soil layer [kg/m^2.s]
Ebare_soil			=	min(Ebare_soil_pot,VavailBare_tm1(1));
Eveg_soil			=	min(Eveg_soil_pot,VavailVeg_tm1(1));
VavailBare_tm1(1)	=	VavailBare_tm1(1) - Ebare_soil;
VavailVeg_tm1(1)	=	VavailVeg_tm1(1) - Eveg_soil;

% Water limitation to ground vegetation transpiration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RfH_Zs_Imp	=	NaN(1,length(Zs)-1);
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,RfH_Zs_ImpL2,~,~,~,~,~,~,~,~,~,~,~]...
	=soil_functions.SoilParametersTotal(Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,Zs(3:end));
RfH_Zs_Imp(1,1:2)	=	[0,0];
RfH_Zs_Imp(1,3:end)	=	RfH_Zs_ImpL2;

% How much water is available per crown area
switch SPARTREE
    case 1	% Tree roots can access all water in the soil (imp, bare, veg)
		Ccrown = [Cveg*fgveg, Cbare*fgbare, Cimp*fgimp, Ctree*4*rad_tree];
		
		Vavail_Veg_tm1_L	=	Ccrown(1)./(Ccrown(1)+Ccrown(1).*Ccrown(4)).*VavailVeg_tm1;
		Vavail_Veg_tm1_H	=	Ccrown(1).*Ccrown(4)./(Ccrown(1)+Ccrown(1).*Ccrown(4)).*VavailVeg_tm1;
		Vavail_Bare_tm1_H	=	Ccrown(2).*Ccrown(4)./(Ccrown(2).*Ccrown(4)).*VavailBare_tm1;
		Vavail_Imp_tm1_H	=	Ccrown(3).*Ccrown(4)./(Ccrown(3).*Ccrown(4)).*VavailImp_tm1;
		
	case 2	% If the tree crown is smaller than the combined vegetated and bare fraction, 
		% then the trees only transpire from these fractions. Otherwise, they
		% also transpire from the impervious ground fraction.
		if (4*rad_tree)<=(fgveg+fgbare)
			Ccrown = [Cveg*fgveg, Cbare*fgbare, Cimp*fgimp, Ctree*4*rad_tree/(fgveg+fgbare)];
		
			Vavail_Veg_tm1_L	=	Ccrown(1)./(Ccrown(1)+Ccrown(1).*Ccrown(4)).*VavailVeg_tm1;
			Vavail_Veg_tm1_H	=	Ccrown(1).*Ccrown(4)./(Ccrown(1)+Ccrown(1).*Ccrown(4)).*VavailVeg_tm1;
			Vavail_Bare_tm1_H	=	Ccrown(2).*Ccrown(4)./(Ccrown(2).*Ccrown(4)).*VavailBare_tm1;
			Vavail_Imp_tm1_H	=	0.*VavailImp_tm1;

		elseif (4*rad_tree)>(fgveg+fgbare)
			Ccrown = [Cveg*fgveg, Cbare*fgbare, Cimp*fgimp, Ctree*1, Ctree*((4*rad_tree)-(fgveg+fgbare))/fgimp];
		
			Vavail_Veg_tm1_L	=	Ccrown(1)./(Ccrown(1)+Ccrown(1).*Ccrown(4)).*VavailVeg_tm1;
			Vavail_Veg_tm1_H	=	Ccrown(1).*Ccrown(4)./(Ccrown(1)+Ccrown(1).*Ccrown(4)).*VavailVeg_tm1;
			Vavail_Bare_tm1_H	=	Ccrown(2).*Ccrown(4)./(Ccrown(2).*Ccrown(4)).*VavailBare_tm1;
			Vavail_Imp_tm1_H	=	Ccrown(3).*Ccrown(5)./(Ccrown(3).*Ccrown(5)).*VavailImp_tm1;
		end
end

% Minimum available and extractable water for plants
Vavail_Veg_tm1_L	= 	sum(min(Vavail_Veg_tm1_L.*RfL_Zs,ExwatVeg_tm1_L*(row/3600/1000)),2);  % [kg/m^2.s]
Vavail_Veg_tm1_H	= 	sum(min(Vavail_Veg_tm1_H.*RfH_Zs,ExwatVeg_tm1_H*(row/3600/1000)),2);  % [kg/m^2.s]
Vavail_Bare_tm1_H	= 	sum(min(Vavail_Bare_tm1_H.*RfH_Zs,ExwatBare_tm1_H*(row/3600/1000)),2);  % [kg/m^2.s]
Vavail_Imp_tm1_H	= 	sum(min(Vavail_Imp_tm1_H.*RfH_Zs_Imp,ExwatImp_tm1_H*(row/3600/1000)),2);  % [kg/m^2.s]

% Water limitation for tree and ground vegetation transpiration
Vavail_plant_tm1_H	=	(Ccrown(1).*Vavail_Veg_tm1_H + Ccrown(2).*Vavail_Bare_tm1_H + Ccrown(3).*Vavail_Imp_tm1_H)./Ccrown(4);
Vavail_plant_tm1_L	=	Vavail_Veg_tm1_L;
TEveg				=	min(TEveg_pot,Vavail_plant_tm1_L);
TEtree				=	min(TEtree_pot,Vavail_plant_tm1_H);

if TEtree_pot>TEtree
	x	=	1;
end

% Evapotranspiration
Ebare		=	Cbare*(Ebare_pond + Ebare_soil);
Eveg		=	Cveg*(Eveg_int+Eveg_pond+Eveg_soil+TEveg);
Etree		=	Ctree*(Etree_int+TEtree);

% Latent heat
LEimp		=	Cimp*(L_heat*Eimp); % [W/m^2]
LEbare_pond	=	Cbare*(L_heat*Ebare_pond);
LEbare_soil	=	Cbare*(L_heat*Ebare_soil);
LEbare		=	Cbare*(LEbare_pond + LEbare_soil);
LEveg_int	=	Cveg*(L_heat*Eveg_int);
LEveg_pond	=	Cveg*(L_heat*Eveg_pond);
LEveg_soil	=	Cveg*(L_heat*Eveg_soil);
LTEveg		=	Cveg*(L_heat*TEveg);
LEveg		=	Cveg*(LEveg_int+LEveg_pond+LEveg_soil+LTEveg);
LEtree_int	=	Ctree*(L_heat*Etree_int);
LTEtree		=	Ctree*(L_heat*TEtree);
LEtree		=	Ctree*(LEtree_int+LTEtree);

