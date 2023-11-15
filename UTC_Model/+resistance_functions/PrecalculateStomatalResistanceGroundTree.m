function[rs_sun_H,rs_shd_H,Ci_sun_H,Ci_shd_H,rs_sun_L,rs_shd_L,Ci_sun_L,Ci_shd_L]=PrecalculateStomatalResistanceGroundTree(...
         TempVec_ittm,Humidity_ittm,ParVegGround,SoilPotW_ittm,CiCO2Leaf_ittm,...
         MeteoData,geometry,FractionsGround,ParTree,...
		 PropOpticalGround,PropOpticalWall,PropOpticalTree,ParVegTree,...
         SunPosition,ViewFactor,ParWindows,BEM_on,...
         rb_L,rb_H,rap_can, rap_Htree_In)


%% Shortwave radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,SWRabs_t,SWRabsDir_t,SWRabsDiff_t,~]...
         =radiation_functions.TotalSWRabsorbed(geometry,FractionsGround,ParTree,...
		 PropOpticalGround,PropOpticalWall,PropOpticalTree,ParVegTree,MeteoData,...
         SunPosition,ViewFactor,ParWindows,BEM_on);
	 
% Tree absorbed: conversion from sphere to horizontal projected area
SWRabs_t.SWRabsTree		=	SWRabs_t.SWRabsTree*4*geometry.radius_tree*pi/(4*geometry.radius_tree);
SWRabsDir_t.SWRabsTree	=	SWRabsDir_t.SWRabsTree*4*geometry.radius_tree*pi/(4*geometry.radius_tree);
SWRabsDiff_t.SWRabsTree	=	SWRabsDiff_t.SWRabsTree*4*geometry.radius_tree*pi/(4*geometry.radius_tree);

SWRdir_abs_tree         = SWRabsDir_t.SWRabsTree;
SWRdiff_abs_tree        = SWRabsDiff_t.SWRabsTree;
SWRdir_abs_groundveg    = SWRabsDir_t.SWRabsGroundVeg;
SWRdiff_abs_groundveg   = SWRabsDiff_t.SWRabsGroundVeg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter specification
Tcanyon			=	TempVec_ittm.TCanyon;
qcanyon			=	Humidity_ittm.CanyonSpecific;
Tveg		    =	TempVec_ittm.TGroundVeg;
Ttree		    =	TempVec_ittm.TTree;
Pre				=	MeteoData.Pre;
Catm_O2			=	MeteoData.Catm_O2;
Catm_CO2		=	MeteoData.Catm_CO2;
fgveg			=	FractionsGround.fveg;
trees			=	ParTree.trees;
LAI_L			=	ParVegGround.LAI;
LAI_H			=	ParVegTree.LAI;
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
Psi_L_tm1		=	SoilPotW_ittm.SoilPotWGroundVeg_L;
Psi_H_tm1		=	SoilPotW_ittm.SoilPotWGroundTot_H;
Ci_sun_H_tm1	=	CiCO2Leaf_ittm.CiCO2LeafTreeSun;
Ci_shd_H_tm1	=	CiCO2Leaf_ittm.CiCO2LeafTreeShd;
Ci_sun_L_tm1	=	CiCO2Leaf_ittm.CiCO2LeafGroundVegSun;
Ci_shd_L_tm1	=	CiCO2Leaf_ittm.CiCO2LeafGroundVegShd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structural parameters
Cveg			=	(fgveg>0);
Ctree			=	(trees==1);

% Vapor pressure saturation and specific humidity at esat
esat_T_canyon	=	611*exp(17.27*(Tcanyon-273.16)/(237.3+(Tcanyon-273.16)));			% vapor pressure saturation at T_canyon [Pa]
e_T_canyon		=	qcanyon*Pre/(0.622+0.378*qcanyon);


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


%% Stomata resistances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Opt_CR	=	optimset('TolFun',1); % [ppm] Numerical tolerance for internal CO2 computation

if Ctree==1 && LAI_H>0
	[rs_sun_H,rs_shd_H,Ci_sun_H,Ci_shd_H,~,~,~,~]=resistance_functions.Canopy_Resistence_An_Evolution(PAR_sun_H,PAR_shd_H,LAI_H,...
		Kopt_H,Knit_H,Fsun_H,Fshd_H,Citm1_sun_H,Citm1_shd_H,...
		Catm_CO2,rap_Htree_In,rb_H,Ttree-273.15,Tcanyon-273.15,Pre/100,Ds_canyon,...
		Psi_H_tm1,Psi_sto_50_H,Psi_sto_00_H,...
		CT_H,Vmax_H,DSE_H,Ha_H,FI_H,Oa,Do_H,a1_H,go_H,e_rel_H,...
		e_relN_H,gmes_H,rjv_H,mSl_H,Sl_H,Opt_CR);
else
	rb_H=Inf;	rs_sun_H=Inf;  rs_shd_H=Inf; Ci_sun_H=0; Ci_shd_H=0;
end

if Cveg==1 && LAI_L>0
	[rs_sun_L,rs_shd_L,Ci_sun_L,Ci_shd_L,~,~,~,~]=resistance_functions.Canopy_Resistence_An_Evolution(PAR_sun_L,PAR_shd_L,LAI_L,...
		Kopt_L,Knit_L,Fsun_L,Fshd_L,Citm1_sun_L,Citm1_shd_L,...
		Catm_CO2,rap_can,rb_L,Tveg-273.15,Tcanyon-273.15,Pre/100,Ds_canyon,...
		Psi_L_tm1,Psi_sto_50_L,Psi_sto_00_L,...
		CT_L,Vmax_L,DSE_L,Ha_L,FI_L,Oa,Do_L,a1_L,go_L,e_rel_L,...
		e_relN_L,gmes_L,rjv_L,mSl_L,Sl_L,Opt_CR);
else
	rb_L=Inf;	rs_sun_L=Inf;  rs_shd_L=Inf; Ci_sun_L=0; Ci_shd_L=0;
end