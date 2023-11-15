function[q_tree_dwn,In_tree,dIn_tree_dt,q_gveg_dwn,In_gveg,dIn_gveg_dt,...
	q_gimp_runoff,In_gimp,dIn_gimp_dt,f_inf_gimp,q_gbare_runoff,In_gbare,dIn_gbare_dt,f_inf_gbare,...
	q_gveg_runoff,In_gveg_pond,dIn_gveg_pond_dt,f_inf_gveg,...
	...
	V_gimp,O_gimp,OS_gimp,Lk_gimp,Psi_s_H_gimp,Psi_s_L_gimp,...
	Exwat_H_gimp,Exwat_L_gimp,Rd_gimp,TEgveg_imp,TEtree_imp,...
	Egimp_soil,dV_dt_gimp,Psi_soil_gimp,Kf_gimp,...
	...
	V_gbare,O_gbare,OS_gbare,Lk_gbare,Psi_s_H_gbare,Psi_s_L_gbare,...
	Exwat_H_gbare,Exwat_L_gbare,Rd_gbare,TEgveg_bare,TEtree_bare,...
	Egbare_Soil,dV_dt_gbare,Psi_soil_gbare,Kf_gbare,...
	...
	V_gveg,O_gveg,OS_gveg,Lk_gveg,Psi_s_H_gveg,Psi_s_L_gveg,...
	Exwat_H_gveg,Exwat_L_gveg,Rd_gveg,TEgveg_veg,TEtree_veg,...
	Egveg_Soil,dV_dt_gveg,Psi_soil_gveg,Kf_gveg,...
	....
	Qin_imp,Qin_bare,Qin_veg,Qin_bare2imp,Qin_bare2veg,Qin_imp2bare,Qin_imp2veg,Qin_veg2imp,Qin_veg2bare,...
	V,O,OS,Lk,Rd,dV_dt,Psi_s_L,Exwat_L,TEgveg_tot,Psi_s_H_tot,Exwat_H,...
	TEtree_tot,EB_TEtree,EB_TEgveg,WBIndv,WBTot,...
	Runoff,Runon_ittm,Etot,DeepGLk,StorageTot]=...
	WaterCanyon(MeteoData,Int_ittm,Owater_ittm,Runon_ittm,Qinlat_ittm,...
	Etree_In,Egveg_In,Egimp_Pond,Egbare_Pond,Egveg_Pond,Egbare_soil,Egveg_soil,TEgveg,TEtree,...
	ParSoilGround,ParInterceptionTree,ParCalculation,ParVegGround,ParVegTree,...
	FractionsGround,geometry,ParTree,Gemeotry_m,Anthropogenic)


%% Assigning input
Rain				=	MeteoData.Rain;
In_gimp_tm1			=	Int_ittm.IntGroundImp;
In_gbare_tm1		=	Int_ittm.IntGroundBare;
In_gveg_tm1			=	Int_ittm.IntGroundVegPlant;
In_gvegpond_tm1		=	Int_ittm.IntGroundVegGround;
In_tree_tm1			=	Int_ittm.IntTree;
Otm1_imp			=	Owater_ittm.OwGroundSoilImp;
Otm1_bare			=	Owater_ittm.OwGroundSoilBare;
Otm1_veg			=	Owater_ittm.OwGroundSoilVeg;
Qin_imp_tm1			=	Qinlat_ittm.Qin_imp;
Qin_bare_tm1		=	Qinlat_ittm.Qin_bare;
Qin_veg_tm1			=	Qinlat_ittm.Qin_veg;
Runon_tm1			=	Runon_ittm.RunonGroundTot;
Wcan				=	Gemeotry_m.Width_canyon;
LAI_g				=	ParVegGround.LAI;
SAI_g				=	ParVegGround.SAI;
LAI_t				=	ParVegTree.LAI;
SAI_t				=	ParVegTree.SAI;
Pcla				=	ParSoilGround.Pcla;
Psan				=	ParSoilGround.Psan;
Porg				=	ParSoilGround.Porg;
Kfc					=	ParSoilGround.Kfc;
Phy					=	ParSoilGround.Phy;
SPAR				=	ParSoilGround.SPAR;
Kbot				=	ParSoilGround.Kbot;
CASE_ROOT_H			=	ParVegTree.CASE_ROOT;
CASE_ROOT_L			=	ParVegGround.CASE_ROOT;
ZR95_H				=	ParVegTree.ZR95;
ZR95_L				=	ParVegGround.ZR95;
ZR50_H				=	ParVegTree.ZR50;
ZR50_L				=	ParVegGround.ZR50;
ZRmax_H				=	ParVegTree.ZRmax;
ZRmax_L				=	ParVegGround.ZRmax;
Zs					=	ParSoilGround.Zs;
In_max_gimp			=	ParSoilGround.In_max_imp;
In_max_gbare		=	ParSoilGround.In_max_bare;
In_max_gvegpond		=	ParSoilGround.In_max_underveg;
Sp_In_tree			=	ParInterceptionTree.Sp_In;
Sp_In_g				=	ParSoilGround.Sp_In;
Kimp				=	ParSoilGround.Kimp;	
dth					=	ParCalculation.dth;
row					=	ParCalculation.row;
Rrootl_H			=	ParVegTree.Rrootl;
Rrootl_L			=	ParVegGround.Rrootl;
PsiL50_H			=	ParVegTree.PsiL50;
PsiL50_L			=	ParVegGround.PsiL50;
PsiX50_H			=	ParVegTree.PsiX50;
PsiX50_L			=	ParVegGround.PsiX50;
SPARTREE			=	ParVegTree.SPARTREE;
Per_runoff			=	FractionsGround.Per_runoff;
fimp				=	FractionsGround.fimp;
fbare				=	FractionsGround.fbare;
fveg				=	FractionsGround.fveg;
r_tree				=	geometry.radius_tree;

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

Cimp			=	(fimp>0);
Cbare			=	(fbare>0);
Cveg			=	(fveg>0);
Ctree			=	(ParTree.trees==1);

%% Water flux compuations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vegetation interception
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trees: water intercepted
[q_tree_dwn,In_tree,dIn_tree_dt,WB_In_tree]=...
	water_functions.Water_Vegetation(Rain,Etree_In,In_tree_tm1,Sp_In_tree,LAI_t,SAI_t,row,dth);

q_tree_dwn	=	q_tree_dwn*Ctree;	% [mm/dth]
In_tree		=	In_tree*Ctree;		% [mm/dth]
dIn_tree_dt	=	dIn_tree_dt*Ctree;	% [mm/dth]
WB_In_tree	=	WB_In_tree*Ctree;	% [mm/dth]

% Water received by any ground fraction including rain and dripping from trees
Rain_ground	=	4*r_tree*Ctree*q_tree_dwn + (1-4*r_tree*Ctree)*Rain; 	% [mm/dth]

% Ground vegetation: interception
% [q_gveg_dwn,In_gveg,dIn_gveg_dt,WB_In_gveg]=...
% 	water_functions.Water_Vegetation(Rain_ground+Anthropogenic.Waterf_canyonVeg,Egveg_In,In_gveg_tm1,Sp_In_g,LAI_g,SAI_g,row,dth);
[q_gveg_dwn,In_gveg,dIn_gveg_dt,WB_In_gveg]=...
	water_functions.Water_Vegetation(Rain_ground,Egveg_In,In_gveg_tm1,Sp_In_g,LAI_g,SAI_g,row,dth);

q_gveg_dwn	=	q_gveg_dwn*Cveg; 	% [mm/dth]
In_gveg		=	In_gveg*Cveg;		% [mm/dth]
dIn_gveg_dt	=	dIn_gveg_dt*Cveg;	% [mm/dth]
WB_In_gveg	=	WB_In_gveg*Cveg;	% [mm/dth]

%% Water ponding on ground 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impervious ground fraction: water ponding
% [~,~,~,f_inf_gimp,~]=...
% 	water_functions.Water_Ground(Rain_ground,Runon_tm1,Egimp_Pond,Otm1_imp(3:end),In_gimp_tm1,In_max_gimp,...
% 	Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
% 	CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L.*NaN,ZR50_H,ZR50_L.*NaN,ZRmax_H,ZRmax_L.*NaN,...
% 	Zs(3:end),dth,row);
%
% f_inf_gimp	=	min(Kimp,f_inf_gimp);

f_inf_gimp	=	Kimp;	% [mm/h]

[q_gimp_runoff,In_gimp,dIn_gimp_dt,f_inf_gimp,WB_In_gimp]=...
	water_functions.Water_Impervious(Rain_ground,Runon_tm1,Egimp_Pond,In_gimp_tm1,dth,row,In_max_gimp,f_inf_gimp);

q_gimp_runoff	=	q_gimp_runoff*Cimp;	% [mm/dth]
In_gimp			=	In_gimp*Cimp;		% [mm/dth]
dIn_gimp_dt		=	dIn_gimp_dt*Cimp;	% [mm/dth]
f_inf_gimp		=	f_inf_gimp*Cimp;	% [mm/h]
WB_In_gimp		=	WB_In_gimp*Cimp;	% [mm/dth]

% Bare ground fraction: water ponding
[q_gbare_runoff,In_gbare,dIn_gbare_dt,f_inf_gbare,WB_In_gbare]=...
	water_functions.Water_Ground(Rain_ground+Anthropogenic.Waterf_canyonBare,Runon_tm1,Egbare_Pond,Otm1_bare,In_gbare_tm1,In_max_gbare,...
	Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L.*NaN,ZR50_H,ZR50_L.*NaN,ZRmax_H,ZRmax_L.*NaN,...
	Zs,dth,row);

q_gbare_runoff	=	q_gbare_runoff*Cbare;	% [mm/dth]
In_gbare		=	In_gbare*Cbare;			% [mm/dth]
dIn_gbare_dt	=	dIn_gbare_dt*Cbare;		% [mm/dth]
f_inf_gbare		=	f_inf_gbare*Cbare;		% [mm/h]
WB_In_gbare		=	WB_In_gbare*Cbare;		% [mm/dth]

% Soil vegetation: water ponding
[q_gveg_runoff,In_gveg_pond,dIn_gveg_pond_dt,f_inf_gveg,WB_Pond_gveg]=...
	water_functions.Water_Ground(q_gveg_dwn+Anthropogenic.Waterf_canyonVeg,Runon_tm1,Egveg_Pond,Otm1_veg,In_gvegpond_tm1,In_max_gvegpond,...
	Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,Zs,dth,row);

q_gveg_runoff	=	q_gveg_runoff*Cveg;		% [mm/dth]
In_gveg_pond	=	In_gveg_pond*Cveg;		% [mm/dth]
dIn_gveg_pond_dt=	dIn_gveg_pond_dt*Cveg;	% [mm/dth]
f_inf_gveg		=	f_inf_gveg*Cveg;		% [mm/h]
WB_Pond_gveg	=	WB_Pond_gveg*Cveg;		% [mm/dth]

%% Water distribution in different soil columns (impervious, bare, vegetated)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transpiration assumptions: the tree can extract water from all soil
% layers. The low vegetation can only take water from the vegetated soil layer
TEtree	=	TEtree;	% Tree evaporation per m^2 tree crown area

switch SPARTREE
    case 1	% Tree roots can access all water in the soil (imp, bare, veg)
		TEtree0_imp		=	TEtree*Ctree*(4*r_tree)*Cimp;
		TEtree0_bare	=	TEtree*Ctree*(4*r_tree)*Cbare;
		TEtree0_veg		=	TEtree*Ctree*(4*r_tree)*Cveg;
	case 2	% If the tree crown is smaller than the combined vegetated and bare fraction, 
		% then the trees only transpire from these fractions. Otherwise, they
		% also transpire from the impervious ground fraction.
		if (4*r_tree)<=(fveg+fbare)
			TEtree0_imp	=	0*Ctree*Cimp;
			TEtree0_bare=	TEtree*(4*r_tree)/(fveg+fbare)*Ctree*Cbare;
			TEtree0_veg	=	TEtree*(4*r_tree)/(fveg+fbare)*Ctree*Cveg;
		elseif (4*r_tree)>(fveg+fbare)
			TEtree0_imp =	((4*r_tree)-(fveg+fbare))*TEtree/fimp*Ctree*Cimp;
			TEtree0_bare=	TEtree*Ctree*Cbare;	
			TEtree0_veg	=	TEtree*Ctree*Cveg;
		end
end

TEgveg		=	TEgveg*Cveg;		% [mm/h]
Egbare_soil	=	Egbare_soil*Cbare;	% [mm/h]
Egveg_soil	=	Egveg_soil*Cveg;	% [mm/h]

% impervious surface
[V_gimp1,O_gimp1,OS_gimp1,Lk_gimp1,Psi_s_H_gimp1,Psi_s_L_gimp1,...
	Exwat_H_gimp1,Exwat_L_gimp1,Rd_gimp1,TEgveg_imp1,TEtree_imp1,...
	Egimp_soil1,dV_dt_gimp1,WB_Soil_gimp1,Psi_Soil_gimp1,Kf_gimp1]=...
	water_functions.Water_Soil(Otm1_imp(3:end),f_inf_gimp,TEtree0_imp,0,0,Qin_imp_tm1(3:end).*0,...
	dth,Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,...
	Rrootl_H,Rrootl_L,PsiL50_H,PsiL50_L,PsiX50_H,PsiX50_L,Zs(3:end)-Zs(3),row);

V_gimp1			=	[NaN, NaN, V_gimp1];
O_gimp1			=	[NaN, NaN, O_gimp1];
Exwat_H_gimp1	=	[NaN, NaN, Exwat_H_gimp1];
Exwat_L_gimp1	=	[NaN, NaN, Exwat_L_gimp1];
Psi_Soil_gimp1	=	[NaN, NaN, Psi_Soil_gimp1];
Kf_gimp1		=	[NaN, NaN, Kf_gimp1];

% bare surface
[V_gbare1,O_gbare1,OS_gbare1,Lk_gbare1,Psi_s_H_gbare1,Psi_s_L_gbare1,Exwat_H_gbare1,Exwat_L_gbare1,Rd_gbare1,TEgveg_bare1,TEtree_bare1,...
	Egbare_Soil1,dV_dt_gbare1,WB_Soil_gbare1,Psi_soil_gbare1,Kf_gbare1]=...
	water_functions.Water_Soil(Otm1_bare,f_inf_gbare,TEtree0_bare,0,Egbare_soil,Qin_bare_tm1.*0,...
	dth,Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,...
	Rrootl_H,Rrootl_L,PsiL50_H,PsiL50_L,PsiX50_H,PsiX50_L,Zs,row);

% vegetated surface
[V_gveg1,O_gveg1,OS_gveg1,Lk_gveg1,Psi_s_H_gveg1,Psi_s_L_gveg1,Exwat_H_gveg1,Exwat_L_gveg1,Rd_gveg1,TEgveg_veg1,TEtree_veg1,...
	Egveg_Soil1,dV_dt_gveg1,WB_Soil_gveg1,Psi_soil_gveg1,Kf_gveg1]=...
	water_functions.Water_Soil(Otm1_veg,f_inf_gveg,TEtree0_veg,TEgveg,Egveg_soil,Qin_veg_tm1.*0,...
	dth,Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,...
	Rrootl_H,Rrootl_L,PsiL50_H,PsiL50_L,PsiX50_H,PsiX50_L,Zs,row);

% LATERAL RICHARDS EQUATION
T_SPAN	=	[0 dth];
OPT_SM	=	odeset('AbsTol',0.05,'MaxStep',dth);
dz		=	diff(Zs); % [mm]  Thickness of the Layers

[~,~,~,Osat,Ohy,nVG,alpVG,Ks_Zs,L,Pe,O33,SPAR,~,~,~,~,~,~,~,~,~,~,~,~,~,~]...
	=soil_functions.SoilParametersTotal(Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,Zs);

V_gimp2		=	zeros(1,length(dz));
V_gbare2	=	zeros(1,length(dz));
V_gveg2		=	zeros(1,length(dz));

% The first two soil layers are only calculated as the exchange between two
% soil columns as the impervious soil column is assumed to be impervious
% there.
for i = 1:2
Vlat1	=	[V_gbare1(i); V_gveg1(i)];

[Tout,Vout]=ode23s(@soil_functions.SOIL_MOISTURES_RICH_COMP_LAT2,T_SPAN,Vlat1,OPT_SM,...
	dz(i),SPAR,Ks_Zs(i),Osat(i),Ohy(i),L(i),Pe(i),O33(i),alpVG(i),nVG(i),...
	Cbare,Cveg,fbare,fveg,Wcan);

V_gimp2(i)	=	NaN;
V_gbare2(i)	=	Vout(end,1);
V_gveg2(i)	=	Vout(end,2);	
end

for i = 3:length(dz)
Vlat1	=	[V_gimp1(i); V_gbare1(i); V_gveg1(i)];

[Tout,Vout]=ode23s(@soil_functions.SOIL_MOISTURES_RICH_COMP_LAT3,T_SPAN,Vlat1,OPT_SM,...
	dz(i),SPAR,Ks_Zs(i),Osat(i),Ohy(i),L(i),Pe(i),O33(i),alpVG(i),nVG(i),...
	Cimp,Cbare,Cveg,fimp,fbare,fveg,Wcan);

V_gimp2(i)	=	Vout(end,1);
V_gbare2(i)	=	Vout(end,2);
V_gveg2(i)	=	Vout(end,3);

%  for jk = 2:length(Tout)
%         %Vlat3	=	[V_gimp2(i); V_gbare2(i); V_gveg2(i)];
%         Vlat3 =  0.5*(Vout(jk,:)+Vout(jk-1,:));
%         [Qin_imp_t(i),Qin_bare_t(i),Qin_veg_t(i),Qin_bare2imp_t(i),Qin_bare2veg_t(i),Qin_imp2bare_t(i),Qin_imp2veg_t(i),Qin_veg2imp_t(i),Qin_veg2bare_t(i)]...
%             =soil_functions.SoilQLateral(Vlat3,dz(i),SPAR,Ks_Zs(i),Osat(i),Ohy(i),L(i),Pe(i),O33(i),alpVG(i),nVG(i),...
%             Cimp,Cbare,Cveg,fimp,fbare,fveg,Wcan);
%         %Qi_out = Qi_out + Qi_out_t*(Tout(jk)-Tout(jk-1)); %%% [mm]
%         Qin_imp(i) = Qin_imp(i) + Qin_imp_t(i)*(Tout(jk)-Tout(jk-1)); %%% [mm]
%         Qin_bare(i) =  Qin_bare(i)  + Qin_bare_t(i)*(Tout(jk)-Tout(jk-1));
%         Qin_veg(i) = Qin_veg(i) +  Qin_veg_t(i)*(Tout(jk)-Tout(jk-1));
% 
% 		Qin_bare2imp(i) = Qin_bare2imp(i) + Qin_bare2imp_t(i)*(Tout(jk)-Tout(jk-1)); %%% [mm]
% 		Qin_bare2veg(i) = Qin_bare2veg(i) + Qin_bare2veg_t(i)*(Tout(jk)-Tout(jk-1)); %%% [mm]
% 		Qin_imp2bare(i) = Qin_imp2bare(i) + Qin_imp2bare_t(i)*(Tout(jk)-Tout(jk-1)); %%% [mm]
% 		Qin_imp2veg(i) = Qin_imp2veg(i) + Qin_imp2veg_t(i)*(Tout(jk)-Tout(jk-1)); %%% [mm]
% 		Qin_veg2imp(i) = Qin_veg2imp(i) + Qin_veg2imp_t(i)*(Tout(jk)-Tout(jk-1)); %%% [mm]
% 		Qin_veg2bare(i) = Qin_veg2bare(i) + Qin_veg2bare_t(i)*(Tout(jk)-Tout(jk-1)); %%% [mm]
% 
%  end

end

% Back compute lateral fluxes
Qin_imp		=	V_gimp2 - V_gimp1;		% [mm/dth]
Qin_bare	=	V_gbare2 - V_gbare1;	% [mm/dth]
Qin_veg		=	V_gveg2 - V_gveg1;		% [mm/dth]

Qin_bare2imp=zeros(1,length(dz)); Qin_bare2veg=zeros(1,length(dz)); Qin_imp2bare=zeros(1,length(dz));
Qin_imp2veg=zeros(1,length(dz)); Qin_veg2imp=zeros(1,length(dz)); Qin_veg2bare=zeros(1,length(dz));

Qin_bare2imp(1:2)	=	[NaN, NaN];
Qin_imp2bare(1:2)	=	[NaN, NaN];
Qin_imp2veg(1:2)	=	[NaN, NaN];
Qin_veg2imp(1:2)	=	[NaN, NaN];

% for i=3:length(dz)
% [Qin_bare2imp(i), Qin_imp2bare(i), Qin_veg2imp(i), Qin_veg2bare(i), Qin_bare2veg(i), Qin_imp2veg(i)]...
% 	=soil_functions.LateralFluxesBackCompute(Qin_imp(i),Qin_bare(i),Qin_veg(i),fimp,fbare,fveg);
% end

% Impervious
[V_gimp2,O_gimp2,OS_gimp2,Psi_Soil_gimp2,Psi_s_H_gimp2,Psi_s_L_gimp2,Exwat_H_gimp2,Exwat_L_gimp2,Kf_gimp2]=...
	soil_functions.SoilMoistureConductivityUpdate(V_gimp2(3:end),Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,Zs(3:end)-Zs(3),...
	Rrootl_H,Rrootl_L,PsiL50_H,PsiL50_L,PsiX50_H,PsiX50_L);

V_gimp2			=	[NaN, NaN, V_gimp2];
O_gimp2			=	[NaN, NaN, O_gimp2];
Exwat_H_gimp2	=	[NaN, NaN, Exwat_H_gimp2];
Exwat_L_gimp2	=	[NaN, NaN, Exwat_L_gimp2];
Psi_Soil_gimp2	=	[NaN, NaN, Psi_Soil_gimp2];
Kf_gimp2		=	[NaN, NaN, Kf_gimp2];

% Bare
[V_gbare2,O_gbare2,OS_gbare2,Psi_soil_gbare2,Psi_s_H_gbare2,Psi_s_L_gbare2,Exwat_H_gbare2,Exwat_L_gbare2,Kf_gbare2]=...
	soil_functions.SoilMoistureConductivityUpdate(V_gbare2,Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,Zs,...
	Rrootl_H,Rrootl_L,PsiL50_H,PsiL50_L,PsiX50_H,PsiX50_L);

% Vegetated
[V_gveg2,O_gveg2,OS_gveg2,Psi_soil_gveg2,Psi_s_H_gveg2,Psi_s_L_gveg2,Exwat_H_gveg2,Exwat_L_gveg2,Kf_gveg2]=...
	soil_functions.SoilMoistureConductivityUpdate(V_gveg2,Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,Zs,...
	Rrootl_H,Rrootl_L,PsiL50_H,PsiL50_L,PsiX50_H,PsiX50_L);

% Change in water volumn
Vtm1_imp	=	(Otm1_imp-Ohy).*dz;  % Water volume in each layer [mm]
Vtm1_imp(1:2)=	[NaN,NaN];
Vtm1_bare	=	(Otm1_bare-Ohy).*dz; % Water volume in each layer [mm]
Vtm1_veg	=	(Otm1_veg-Ohy).*dz;  % Water volume in each layer [mm]

dV_dt_gimpTot		=	nansum(V_gimp2,2)-nansum(Vtm1_imp,2);	% [mm/dth]
dV_dt_gbareTot		=	nansum(V_gbare2,2)-nansum(Vtm1_bare,2);	% [mm/dth]
dV_dt_gvegTot		=	nansum(V_gveg2,2)-nansum(Vtm1_veg,2);	% [mm/dth]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface water balance: Tree, impervious, bare, vegetated ground
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WBsurf_tree	=	Ctree*Rain - Etree_In*dth*3600*1000/row - q_tree_dwn - dIn_tree_dt; % Is not used [mm/dth]

WBsurf_imp	=	Cimp*Rain_ground + Cimp*Runon_tm1 - Egimp_Pond*dth*3600*1000/row - q_gimp_runoff - f_inf_gimp*dth - dIn_gimp_dt;	% [mm/dth]

WBsurf_bare	=	Cbare*Rain_ground + Cbare*Runon_tm1 + Anthropogenic.Waterf_canyonBare - Egbare_Pond*dth*3600*1000/row - ...
				f_inf_gbare*dth  - q_gbare_runoff - dIn_gbare_dt;	% [mm/dth]
					
WBsurf_veg	=	Cveg*Rain_ground + Cveg*Runon_tm1 + Anthropogenic.Waterf_canyonVeg...
				- (Egveg_In + Egveg_Pond)*dth*3600*1000/row - ...
				f_inf_gveg*dth - q_gveg_runoff - dIn_gveg_dt - dIn_gveg_pond_dt;	% [mm/dth]
					
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soil water balance: impervious, bare, vegetated ground
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WBsoil_imp	=	f_inf_gimp*dth + nansum(Qin_imp) - (TEgveg_imp1 + TEtree_imp1 + Egimp_soil1)*dth*3600*1000/row...
				- Lk_gimp1*dth - Rd_gimp1 - dV_dt_gimpTot;	% [mm/dth]
						
WBsoil_bare	=	f_inf_gbare*dth + nansum(Qin_bare) - (TEgveg_bare1 + TEtree_bare1 + Egbare_Soil1)*dth*3600*1000/row...
				- Lk_gbare1*dth - Rd_gbare1 - dV_dt_gbareTot;	% [mm/dth]
						
WBsoil_veg	=	f_inf_gveg*dth + nansum(Qin_veg) - (TEgveg_veg1 + TEtree_veg1 + Egveg_Soil1)*dth*3600*1000/row...
				- Lk_gveg1*dth - Rd_gveg1 - dV_dt_gvegTot;	% [mm/dth]
			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ground water balance: impervious, bare, vegetated ground
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%						
WBimp_tot	=	Cimp*Rain_ground + Cimp*Runon_tm1 + nansum(Qin_imp)...
				- (Egimp_Pond+TEgveg_imp1 + TEtree_imp1 + Egimp_soil1)*dth*3600*1000/row...
				- q_gimp_runoff - dIn_gimp_dt - Lk_gimp1*dth - Rd_gimp1 - dV_dt_gimpTot; 	% [mm/dth]

WBbare_tot	=	Cbare*Rain_ground + Cbare*Runon_tm1 + nansum(Qin_bare)+ Anthropogenic.Waterf_canyonBare...
				- (Egbare_Pond + TEgveg_bare1 + TEtree_bare1 + Egbare_Soil1)*dth*3600*1000/row - ...
				- q_gbare_runoff - dIn_gbare_dt - Lk_gbare1*dth - Rd_gbare1 - dV_dt_gbareTot; 	% [mm/dth]
									
WBveg_tot	=	Cveg*Rain_ground + Cveg*Runon_tm1 + nansum(Qin_veg) + Anthropogenic.Waterf_canyonVeg...
				- (Egveg_In + Egveg_Pond + TEgveg_veg1 + TEtree_veg1 + Egveg_Soil1)*dth*3600*1000/row - ...
				- q_gveg_runoff - dIn_gveg_dt - dIn_gveg_pond_dt - Lk_gveg1*dth - Rd_gveg1 - dV_dt_gvegTot;	 	% [mm/dth]				
		
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%				
Runoff		=	Per_runoff*(fimp*(q_gimp_runoff + Rd_gimp1) + fbare*(q_gbare_runoff + Rd_gbare1)...
				+ fveg*(q_gveg_runoff + Rd_gveg1)); 	% [mm/dth]
					
Runon_ittm		=	(1 - Per_runoff)*(fimp*(q_gimp_runoff + Rd_gimp1) + fbare*(q_gbare_runoff + Rd_gbare1)...
				+ fveg*(q_gveg_runoff + Rd_gveg1)); 	% [mm/dth]

Etot		=	(fimp*(Egimp_Pond + TEgveg_imp1 + TEtree_imp1 + Egimp_soil1)...
				+ fbare*(Egbare_Pond + Egbare_Soil1 + TEgveg_bare1 + TEtree_bare1)...
				+ fveg*(Egveg_In + Egveg_Pond + Egveg_Soil1 + TEgveg_veg1 + TEtree_veg1)...
				+ Ctree*(4*r_tree)*(Etree_In))*dth*3600*1000/row; 	% [mm/dth]

DeepGLk		=	(fimp*Lk_gimp1 + fbare*Lk_gbare1 +fveg*Lk_gveg1)*dth; 	% [mm/dth]

StorageTot	=	fimp*(dIn_gimp_dt+dV_dt_gimpTot) + fbare*(dIn_gbare_dt+dV_dt_gbareTot) + ...
				fveg*(dIn_gveg_dt + dIn_gveg_pond_dt + dV_dt_gvegTot) + ...
				(4*r_tree)*dIn_tree_dt; 	% [mm/dth]

WBcanyon_flux	=	Rain + Runon_tm1 + fveg*Anthropogenic.Waterf_canyonVeg + fbare*Anthropogenic.Waterf_canyonBare...
					- Etot - DeepGLk - Runoff - Runon_ittm - StorageTot; %[mm]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
WBtree_level	=	(4*r_tree)*(Ctree*Rain - Etree_In*dth*3600*1000/row - q_tree_dwn - dIn_tree_dt); % Is not used 	% [mm/dth]

WBground_level	=	fimp*(Rain_ground + Runon_tm1-dIn_gimp_dt - Egimp_Pond*dth*3600*1000/row - q_gimp_runoff - Rd_gimp1 - f_inf_gimp*dth)...
					+ fbare*(Rain_ground + Runon_tm1 + Anthropogenic.Waterf_canyonBare -dIn_gbare_dt - Egbare_Pond*dth*3600*1000/row - q_gbare_runoff - Rd_gbare1 - f_inf_gbare*dth)...
					+ fveg*(Rain_ground + Runon_tm1 + Anthropogenic.Waterf_canyonVeg -dIn_gveg_dt - dIn_gveg_pond_dt -(Egveg_In + Egveg_Pond)*dth*3600*1000/row - q_gveg_runoff - Rd_gveg1 - f_inf_gveg*dth);	% [mm/dth]

WBsoil_level	=	fimp*(f_inf_gimp*dth + 3600*1000/row*dth*(- TEgveg_imp1 - TEtree_imp1 - Egimp_soil1) - Lk_gimp1*dth - dV_dt_gimpTot)...
					+ fbare*(f_inf_gbare*dth + 3600*1000/row*dth*(- TEgveg_bare1 - TEtree_bare1 - Egbare_Soil1) - Lk_gbare1*dth - dV_dt_gbareTot)...
					+ fveg*(f_inf_gveg*dth + 3600*1000/row*dth*(- TEgveg_veg1 - TEtree_veg1 - Egveg_Soil1) - Lk_gveg1*dth - dV_dt_gvegTot);	% [mm/dth]

					
WBcanyon_level	=	Rain + Runon_tm1....
					+ (4*r_tree)*(- Etree_In*dth*3600*1000/row - dIn_tree_dt)...
					+ fimp*(-dIn_gimp_dt - Egimp_Pond*dth*3600*1000/row - q_gimp_runoff - Rd_gimp1 - 3600*1000*dth/row*(TEgveg_imp1 + TEtree_imp1 + Egimp_soil1) - Lk_gimp1*dth - dV_dt_gimpTot)...
					+ fbare*(Anthropogenic.Waterf_canyonBare-dIn_gbare_dt - Egbare_Pond*dth*3600*1000/row - q_gbare_runoff - Rd_gbare1 - 3600*1000/row*dth*(TEgveg_bare1 + TEtree_bare1 + Egbare_Soil1) - Lk_gbare1*dth - dV_dt_gbareTot)...
					+ fveg*(Anthropogenic.Waterf_canyonVeg -dIn_gveg_dt - dIn_gveg_pond_dt -(Egveg_In + Egveg_Pond)*dth*3600*1000/row - q_gveg_runoff - Rd_gveg1 - 3600*1000/row*dth*(TEgveg_veg1 + TEtree_veg1 + Egveg_Soil1) - Lk_gveg1*dth - dV_dt_gvegTot);
					% [mm/dth]
					
				
% Assigning the write variables
V_gimp					=	V_gimp2;
O_gimp					=	O_gimp2;
V_gimp(isnan(V_gimp))	=	0;
O_gimp(isnan(O_gimp))	=	0;%Ohy(isnan(O_gimp));
OS_gimp			=	OS_gimp2;
Lk_gimp			=	Lk_gimp1;
Psi_s_H_gimp	=	Psi_s_H_gimp2;
Psi_s_L_gimp	=	Psi_s_L_gimp2;
Exwat_H_gimp	=	Exwat_H_gimp2;
Exwat_L_gimp	=	Exwat_L_gimp2;
Exwat_H_gimp(isnan(Exwat_H_gimp))	=	0;
Exwat_L_gimp(isnan(Exwat_L_gimp))	=	0;
Rd_gimp			=	Rd_gimp1;
TEgveg_imp		=	TEgveg_imp1;
Egimp_soil		=	Egimp_soil1;
dV_dt_gimp		=	dV_dt_gimpTot;
Psi_soil_gimp	=	Psi_Soil_gimp2;
Kf_gimp			=	Kf_gimp2;

V_gbare			=	V_gbare2;
O_gbare			=	O_gbare2;
OS_gbare		=	OS_gbare2;
Lk_gbare		=	Lk_gbare1;
Psi_s_H_gbare	=	Psi_s_H_gbare2;
Psi_s_L_gbare	=	Psi_s_L_gbare2;
Exwat_H_gbare	=	Exwat_H_gbare2;
Exwat_L_gbare	=	Exwat_L_gbare2;
Rd_gbare		=	Rd_gbare1;
TEgveg_bare		=	TEgveg_bare1;
Egbare_Soil		=	Egbare_Soil1;
dV_dt_gbare		=	dV_dt_gbareTot;
Psi_soil_gbare	=	Psi_soil_gbare2;
Kf_gbare		=	Kf_gbare2;

V_gveg			=	V_gveg2;
O_gveg			=	O_gveg2;
OS_gveg			=	OS_gveg2;
Lk_gveg			=	Lk_gveg1;
Psi_s_H_gveg	=	Psi_s_H_gveg2;
Psi_s_L_gveg	=	Psi_s_L_gveg2;
Exwat_H_gveg	=	Exwat_H_gveg2;
Exwat_L_gveg	=	Exwat_L_gveg2;
Rd_gveg			=	Rd_gveg1;
TEgveg_veg		=	TEgveg_veg1;
Egveg_Soil		=	Egveg_Soil1;
dV_dt_gveg		=	dV_dt_gvegTot;
Psi_soil_gveg	=	Psi_soil_gveg2;
Kf_gveg			=	Kf_gveg2;

WB_Soil_gimp	=	WB_Soil_gimp1;
WB_Soil_gbare	=	WB_Soil_gbare1;
WB_Soil_gveg	=	WB_Soil_gveg1;

				
				
% Average properties and rescaling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V	=	nansum([fimp.*V_gimp;fbare.*V_gbare;fveg.*V_gveg],1);
O	=	nansum([fimp.*O_gimp;fbare.*O_gbare;fveg.*O_gveg],1);
OS	=	nansum([fimp.*OS_gimp;fbare.*OS_gbare;fveg.*OS_gveg],1);
Lk	=	nansum([fimp.*Lk_gimp;fbare.*Lk_gbare;fveg.*Lk_gveg],1);
Rd	=	nansum([fimp.*Rd_gimp;fbare.*Rd_gbare;fveg.*Rd_gveg],1);
dV_dt=	nansum([fimp.*dV_dt_gimp;fbare.*dV_dt_gbare;fveg.*dV_dt_gveg],1);

V(1,1:2)=	V(1,1:2)./(fbare+fveg);
O(1,1:2)=	O(1,1:2)./(fbare+fveg);
V(isnan(V))	=	0;
O(isnan(O))	=	0;%Ohy(isnan(O));

Psi_s_L		=	Psi_s_L_gveg;
Exwat_L		=	Exwat_L_gveg;
TEgveg_tot	=	TEgveg_veg;

Psi_s_H_tot	=	Psi_s_H_gveg;
Exwat_H		=	nansum([fimp.*Exwat_H_gimp;fbare.*Exwat_H_gbare;fveg.*Exwat_H_gveg],1);

switch SPARTREE
    case 1	% Tree roots can access all water in the soil (imp, bare, veg)
		TEtree_imp	=	TEtree_imp1/(4*r_tree)*Cimp;
		TEtree_bare	=	TEtree_bare1/(4*r_tree)*Cbare;
		TEtree_veg	=	TEtree_veg1/(4*r_tree)*Cveg;
	case 2	% If the tree crown is smaller than the combined vegetated and bare fraction, 
		% then the trees only transpire from these fractions. Otherwise, they
		% also transpire from the impervious ground fraction.
		if (4*r_tree)<=(fveg+fbare)
			TEtree_imp	=	0*Cimp;
			TEtree_bare	=	TEtree_bare1/(4*r_tree)*Cbare;
			TEtree_veg	=	TEtree_veg1/(4*r_tree)*Cveg;
		elseif (4*r_tree)>(fveg+fbare)
			TEtree_imp	=	TEtree_imp1/((4*r_tree)-(fveg+fbare))*fimp*Cimp;
			TEtree_bare	=	TEtree_bare1*Cbare;	
			TEtree_veg	=	TEtree_veg1*Cveg;
		end
end

TEtree_tot	=	(fimp*TEtree_imp + fbare*TEtree_bare + fveg*TEtree_veg);

EB_TEtree	=	TEtree_tot-TEtree;
EB_TEgveg	=	TEgveg_tot-TEgveg;

% Water balance structs	
WBIndv	=	struct('WB_In_tree',WB_In_tree,'WB_In_gveg',WB_In_gveg,...
			'WB_In_gimp',WB_In_gimp,'WB_In_gbare',WB_In_gbare,...
			'WB_Pond_gveg',WB_Pond_gveg,'WB_Soil_gimp',WB_Soil_gimp,...
			'WB_Soil_gbare',WB_Soil_gbare,'WB_Soil_gveg',WB_Soil_gveg);

WBTot	=	struct('WBsurf_tree',WBsurf_tree,'WBsurf_imp',WBsurf_imp,...
			'WBsurf_bare',WBsurf_bare,'WBsurf_veg',WBsurf_veg,...
			'WBsoil_imp',WBsoil_imp,'WBsoil_bare',WBsoil_bare,...
			'WBsoil_veg',WBsoil_veg,'WBimp_tot',WBimp_tot,...
			'WBbare_tot',WBbare_tot,'WBveg_tot',WBveg_tot,...
			'WBcanyon_flux',WBcanyon_flux,...
			'WBtree_level',WBtree_level,'WBground_level',WBground_level,...
			'WBsoil_level',WBsoil_level,'WBcanyon_level',WBcanyon_level);

end





