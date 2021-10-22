function[Gemeotry_m,ParTree,geometry,FractionsRoof,FractionsGround,...
	WallLayers,ParSoilRoof,ParSoilGround,ParInterceptionTree,...
	PropOpticalRoof,PropOpticalGround,PropOpticalWall,PropOpticalTree,...
	ParThermalRoof,ParThermalGround,ParThermalWall,ParThermalTree,...
	ParVegRoof,ParVegGround,ParVegTree,Person]=Data_UEHM_site(MeteoData,varargin)

% Assign vlaues in case there is a varying LAI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
doy		=	day(MeteoData.Time,'dayofyear'); % Day of year

if isstruct(varargin{1, 2})
    LAIvarying  = 1;
    LAI_Rvar	=   varargin{1, 2}.LAI_R(doy);
    LAI_Gvar	=   varargin{1, 2}.LAI_G(doy);
    LAI_Tvar	=	varargin{1, 2}.LAI_T(doy);
    
    LAI_Rvar(LAI_Rvar==0) = 10^-6; % Somehow the solver fails if LAI_T is zero. But it can be very low
    LAI_Gvar(LAI_Gvar==0) = 10^-6;
    LAI_Tvar(LAI_Tvar==0) = 10^-6;
else 
    LAIvarying = 0;
    LAI_Rvar = 0; LAI_Gvar = 0; LAI_Tvar = 0;
end
    

%% GEOMETRY OF URBAN AREA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Height_canyon	=	6.5;				% Canyon (building) height (m)
Width_canyon	=	5.78;				% Canyon (road) width (m)
Width_roof		=	4.73;				% Roof width (m), calculated from the land cover fraction and the street width.
Radius_tree		=	0.289;				% Tree-crown radius (m), Calculated out of rescaled tree fraction and street width (assuming two uniform strips of tree rows)
Height_tree		=	5-Radius_tree;		% Tree height (m)
Distance_tree	=	0.5+Radius_tree;	% Tree-to-wall distance (m)

Hcan_max	=	NaN;	% Maximum height of roughness elements (buidlings), (m)
Hcan_std	=	NaN;	% Standard deviation of roughness elements (buildings), (m)

trees   =	1;		% Easy switch to include (=1) and exclude (=0) trees in the urban canyon
ftree	=	1;		% DO NOT CHANGE: Tree fraction along canyon axis

if isnan(Radius_tree)==1	% Tree radius cannot be NaN
	Radius_tree	=	0;
end

hcanyon			=	Height_canyon/Width_canyon;		% normalized canyon height(-)
wcanyon			=	Width_canyon/Width_canyon;		% normalized canyon width (-)
wroof			=	Width_roof/Width_canyon;		% normalized roof width (-)
htree			=	Height_tree/Width_canyon;		% normalized tree height (-)
radius_tree		=	Radius_tree/Width_canyon;		% normalized tree radius (-)
distance_tree	=	Distance_tree/Width_canyon;		% normalized tree-to-wall distance (-)
ratio			=	hcanyon/wcanyon;				% height to width ratio (-)

wcanyon_norm	=	wcanyon/(wcanyon+wroof);		% normalized canyon width overall (-)
wroof_norm		=	wroof/(wcanyon+wroof);			% normalized roof width overall (-)

Gemeotry_m	=	struct('Height_canyon',Height_canyon,'Width_canyon',Width_canyon,...
				'Width_roof',Width_roof,'Height_tree',Height_tree,...
				'Radius_tree',Radius_tree,'Distance_tree',Distance_tree,...
				'Hcan_max',Hcan_max,'Hcan_std',Hcan_std);

ParTree		=	struct('trees',trees,'ftree',ftree);

geometry	=	struct('hcanyon',hcanyon,'wcanyon',wcanyon,'wroof',wroof,...
				'htree',htree,'radius_tree',radius_tree,'distance_tree',distance_tree,...
				'ratio',ratio,'wcanyon_norm',wcanyon_norm,'wroof_norm',wroof_norm);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SURFACE FRACTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fveg_R	=	0.5;	% Vegetated roof raction (-)
fimp_R	=	0.5;	% Impervious roof fraction (-)

Per_runoff_R	=	1;		% Percentage of excess water that leaves the system as runoff, needs to be between 0-1 [-]

FractionsRoof	=	struct('fveg',fveg_R,'fimp',fimp_R,'Per_runoff',Per_runoff_R);

% GROUND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fveg_G	=	0.545;	% Vegetated ground raction (-)
fbare_G	=	0.0;	% Bare ground raction (-)
fimp_G	=	0.455;	% Impervious ground raction (-)

Per_runoff_G	=	0.9;	% Percentage of excess water that leaves the system as runoff, needs to be between 0-1 [-]

FractionsGround	=	struct('fveg',fveg_G,'fbare',fbare_G,'fimp',fimp_G,'Per_runoff',Per_runoff_G);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROOF VEGETATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General
LAI_RC      =   2.5;
LAI_R		=	LAIvarying*LAI_Rvar + (1-LAIvarying) * LAI_RC;         % Leaf area index for the roof vegetation (-)
SAI_R		=	0.001;		% Stem area index for the roof vegetation (-)
hc_R		=	0.15;		% canopy height roof vegetation	(m)
h_disp_R	=	2/3*hc_R;	% Zero plane displacement height of roof vegetation [m]
d_leaf_R	=	0.8;		% Leaf dimension of roof vegetation [cm]

% Roof water uptake
CASE_ROOT_R	=	1;		% Type of Root Profile
ZR95_R		=	70;		% Root depth 95 percentile [mm]
ZR50_R		=	NaN;	% Root depth 50 percentile [mm]
ZRmax_R		=	NaN;	% Maximum Root depth [mm]
Rrootl_R	=	3800;	% Root length index [m root / m^2 PFT]
PsiL50_R	=	-4.0;	% Water Potential at 50% loss conductivity [MPa]
PsiX50_R	=	-4.50;	% Water potential at 50 of xylem hydraulic conductivity and limit for water extraction from soil [MPa]

% Photosynthesis and Transpiration
FI_R		=	0.081;	% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_R		=	1000;	% [Pa] Empirical coefficient for the role of vapor pressure in the biochemical model of photosynthesis
a1_R		=	6;		% [-] Empirical parameter connecting stomatal aperture and net assimilaton
go_R		=	0.01;	% [mol / s m^2] minimal Stomatal Conductance
CT_R		=	3;		% --> 'CT' == 3  'CT' ==  4 Photosyntesis Typology for Plants, Photosynthetic pathway C3 or C4
DSE_R		=	0.656;	% [kJ/mol] Activation Energy - Plant Dependent, Activation Energy in Photosynthesis for Rubisco Capacity
Ha_R		=	55;		% [kJ / mol K]  entropy factor - Plant Dependent, Activation energy.
gmes_R		=	Inf;	% [mol CO2/s m2] Mesophyll conductance, not used
rjv_R		=	2.4;	% [?mol Eq/ ?mol CO2] Scaling factor between Jmax and Vmax
Kopt_R		=	0.5;	% [-] optical depth of direct beam perunit plant area
Knit_R		=	0.15;	% [-] Canopy nitrogen decay coefficient
Vmax_R		=	68;		% [?mol CO2/ m2 s] Maximum Rubisco capacity at 25°C leaf level
mSl_R		=	0.0;
e_rel_R		=	1;		% [-] Relative Efficiency of the photosynthesis apparatus due to Age/Day-length
e_relN_R	=	1;		% [-] Relative efficiency of the photosynthesis apparatus due to N limitations
Psi_sto_00_R=	-0.5;	% [MPa]  Water Potential at PLCs loss conductivity
Psi_sto_50_R=	-3.0;	% [MPa]  Water Potential at 50% loss conductivity
Sl_R		=	0.035;	% specific leaf area of  biomass [m^2 /gC]

ParVegRoof	=	struct('LAI',LAI_R,'SAI',SAI_R,'hc',hc_R,'h_disp',h_disp_R,...
				'd_leaf',d_leaf_R,'CASE_ROOT',CASE_ROOT_R,'ZR95',ZR95_R,'ZR50',ZR50_R,...
				'ZRmax',ZRmax_R,'Rrootl',Rrootl_R,'PsiL50',PsiL50_R,'PsiX50',PsiX50_R,...
				'FI',FI_R,'Do',Do_R,'a1',a1_R,'go',go_R,'CT',CT_R,'DSE',DSE_R,'Ha',Ha_R,...
				'gmes',gmes_R,'rjv',rjv_R,'Kopt',Kopt_R,'Knit',Knit_R,'Vmax',Vmax_R,...
				'mSl',mSl_R,'e_rel',e_rel_R,'e_relN',e_relN_R,'Psi_sto_00',Psi_sto_00_R,...
				'Psi_sto_50',Psi_sto_50_R,'Sl',Sl_R);
				
% GROUND VEGETATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRASS
% General
LAI_GC      =   2.5;
LAI_G		=	LAIvarying*LAI_Gvar + (1-LAIvarying) * LAI_GC; 		% Leaf area index for the roof vegetation (-)
SAI_G		=	0.001;		% Stem area index for the roof vegetation (-)
hc_G		=	0.15;		% canopy height roof vegetation (m)
h_disp_G	=	2/3*hc_G;	% Zero plane displacement height of roof vegetation [m]
d_leaf_G	=	0.8;		% Leaf dimension of roof vegetation [cm]

% Roof water uptake
CASE_ROOT_G	=	1;		% Type of Root Profile
ZR95_G		=	250;	% Root depth 95 percentile [mm]
ZR50_G		=	NaN;	% Root depth 50 percentile [mm]
ZRmax_G		=	NaN;	% Maximum Root depth [mm]
Rrootl_G	=	3800;	% Root length index [m root / m^2 PFT]
PsiL50_G	=	-4.0;	% [MPa]  Water Potential at 50% loss conductivity
PsiX50_G	=	-4.50;	% Water potential at 50 of xylem hydraulic conductivity and limit for water extraction from soil [MPa]

% Photosynthesis and Transpiration
FI_G		=	0.081;	% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_G		=	1000;	% [Pa] Empirical coefficient for the role of vapor pressure in the biochemical model of photosynthesis
a1_G		=	6;		% [-] Empirical parameter connecting stomatal aperture and net assimilaton
go_G		=	0.01;	% [mol / s m^2] minimal Stomatal Conductance
CT_G		=	3;		% --> 'CT' == 3  'CT' ==  4 Photosyntesis Typology for Plants, Photosynthetic pathway C3 or C4
DSE_G		=	0.656;	% [kJ/mol] Activation Energy - Plant Dependent, Activation Energy in Photosynthesis for Rubisco Capacity
Ha_G		=	55;		% [kJ / mol K]  entropy factor - Plant Dependent, Activation energy.
gmes_G		=	Inf;	% [mol CO2/s m2] Mesophyll conductance, not used
rjv_G		=	2.4;	% [?mol Eq/ ?mol CO2] Scaling factor between Jmax and Vmax
Kopt_G		=	0.5;	% [-] optical depth of direct beam perunit plant area
Knit_G		=	0.15;	% [-] Canopy nitrogen decay coefficient
Vmax_G		=	68;		% [?mol CO2/ m2 s] Maximum Rubisco capacity at 25°C leaf level
mSl_G		=	0.0;
e_rel_G		=	1;		% [-] Relative Efficiency of the photosynthesis apparatus due to Age/Day-length
e_relN_G	=	1;		% [-] Relative efficiency of the photosynthesis apparatus due to N limitations
Psi_sto_00_G=	-0.5;	% [MPa]  Water Potential at PLCs loss conductivity
Psi_sto_50_G=	-3.0;	% [MPa]  Water Potential at 50% loss conductivity
Sl_G		=	0.035;	% 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]

ParVegGround	=	struct('LAI',LAI_G,'SAI',SAI_G,'hc',hc_G,'h_disp',h_disp_G,...
					'd_leaf',d_leaf_G,'CASE_ROOT',CASE_ROOT_G,'ZR95',ZR95_G,'ZR50',ZR50_G,...
					'ZRmax',ZRmax_G,'Rrootl',Rrootl_G,'PsiL50',PsiL50_G,'PsiX50',PsiX50_G,...
					'FI',FI_G,'Do',Do_G,'a1',a1_G,'go',go_G,'CT',CT_G,'DSE',DSE_G,'Ha',Ha_G,...
					'gmes',gmes_G,'rjv',rjv_G,'Kopt',Kopt_G,'Knit',Knit_G,'Vmax',Vmax_G,...
					'mSl',mSl_G,'e_rel',e_rel_G,'e_relN',e_relN_G,'Psi_sto_00',Psi_sto_00_G,...
					'Psi_sto_50',Psi_sto_50_G,'Sl',Sl_G);
					
% TREE VEGETATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General
LAI_TC      =   5;
LAI_T		=	LAIvarying*LAI_Tvar + (1-LAIvarying) * LAI_TC; 		% Leaf area index for the ground vegetation (-)
SAI_T		=	0.2;	% Stem area index for the ground vegetation (-)
d_leaf_T	=	4;		% Leaf dimension of ground vegetation [cm]

% Tree water uptake
CASE_ROOT_T	=	1;		% Type of Root Profile
ZR95_T		=	1000;	% Root depth 95 percentile [mm]
ZR50_T		=	NaN;	% Root depth 50 percentile [mm]
ZRmax_T		=	NaN;	% Maximum Root depth [mm]
Rrootl_T	=	4000;	% Root length index [m root / m^2 PFT]
PsiL50_T	=	-3.0;	% [MPa]  Water Potential at 50% loss conductivity
PsiX50_T	=	-4.5;	% Water potential at 50 of xylem hydraulic conductivity and limit for water extraction from soil

% Photosynthesis and Transpiration
FI_T		=	0.081;	% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_T		=	1000;	% [Pa] Empirical coefficient for the role of vapor pressure in the biochemical model of photosynthesis
a1_T		=	9;		% [-] Empirical parameter connecting stomatal aperture and net assimilaton
go_T		=	0.01;	% [mol / s m^2] minimal Stomatal Conductance
CT_T		=	3;		% --> 'CT' == 3  'CT' ==  4 Photosyntesis Typology for Plants, Photosynthetic pathway C3 or C4
DSE_T		=	0.649;	% [kJ/mol] Activation Energy - Plant Dependent, Activation Energy in Photosynthesis for Rubisco Capacity
Ha_T		=	76;		% [kJ / mol K]  entropy factor - Plant Dependent, Activation energy.
gmes_T		=	Inf;	% [mol CO2/s m2] Mesophyll conductance, not used
rjv_T		=	2.4;	% [?mol Eq/ ?mol CO2] Scaling factor between Jmax and Vmax
Kopt_T		=	0.5;	% [-] optical depth of direct beam perunit plant area ???
Knit_T		=	0.35;	% [-] Canopy nitrogen decay coefficient
Vmax_T		=	66;		% [?mol CO2/ m2 s] Maximum Rubisco capacity at 25°C leaf level
mSl_T		=	0.0;
e_rel_T		=	1;		% [-] Relative Efficiency of the photosynthesis apparatus due to Age/Day-length
e_relN_T	=	1;		% [-] Relative efficiency of the photosynthesis apparatus due to N limitations
Psi_sto_00_T=	-0.5;	% [MPa]  Water Potential at PLCs loss conductivity
Psi_sto_50_T=	-2.2;	% [MPa]  Water Potential at 50% loss conductivity
Sl_T		=	0.024;	% 0.05 -0.005 [m^2 gC] specific leaf area of  biomass [m^2 /gC]
SPARTREE	=	2;		% Tree root distribution: 1 = Tree roots can access all water in the soil (imp, bare, veg) equally
						% 2 =  If the tree crown is smaller than the combined vegetated and bare fraction, 
						% then the trees only transpire from these fractions. Otherwise, they
						% also transpire from the impervious ground fraction.

ParVegTree	=	struct('LAI',LAI_T,'SAI',SAI_T,'d_leaf',d_leaf_T,'CASE_ROOT',CASE_ROOT_T,...
				'ZR95',ZR95_T,'ZR50',ZR50_T,'ZRmax',ZRmax_T,'Rrootl',Rrootl_T,...
				'PsiL50',PsiL50_T,'PsiX50',PsiX50_T,'FI',FI_T,'Do',Do_T,'a1',a1_T,...
				'go',go_T,'CT',CT_T,'DSE',DSE_T,'Ha',Ha_T,'gmes',gmes_T,'rjv',rjv_T,...
				'Kopt',Kopt_T,'Knit',Knit_T,'Vmax',Vmax_T,'mSl',mSl_T,'e_rel',e_rel_T,...
				'e_relN',e_relN_T,'Psi_sto_00',Psi_sto_00_T,'Psi_sto_50',Psi_sto_50_T,...
				'Sl',Sl_T,'SPARTREE',SPARTREE);
            
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OPTICAL PROPERTIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROOF OPTICAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aveg_R		=	0.2;	% Roof vegetation surface albedo (-)
aimp_R		=	0.15;	% Roof impervious albedo (-)
albedo_R	=	fveg_R*aveg_R+fimp_R*aimp_R;	% equivalent roof surface albedo (-)

eveg_R		=	1 - exp(-(LAI_R+SAI_R));	% Roof vegetation emissivity (-) 
eimp_R		=	0.95;	% Roof impervious emissivity (-)
emissivity_R=	fveg_R*eveg_R+fimp_R*eimp_R;	% equivalent roof surface emissivity (-)

PropOpticalRoof	=	struct('aveg',aveg_R,'aimp',aimp_R,'albedo',albedo_R,...
					'eveg',eveg_R,'eimp',eimp_R,'emissivity',emissivity_R);

% GROUND OPTICAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aveg_G		=	0.2;	% Ground vegetation surface albedo (-)
abare_G		=	0.15;	% Ground vegetation surface albedo (-)
aimp_G		=	0.1;	% Ground impervious albedo (-)
albedo_G	=	fveg_G*aveg_G + fbare_G*abare_G + fimp_G*aimp_G;	% equivalent ground surface albedo (-)

eveg_G		=	1 - exp(-(LAI_G+SAI_G));	% Ground vegetation emissivity (-) 
ebare_G		=	0.95;	% Ground vegetation surface emissivity (-)
eimp_G		=	0.95;	% Ground impervious emissivity (-)
emissivity_G=	fveg_G*eveg_G + fbare_G*ebare_G + fimp_G*eimp_G;	% equivalent ground surface emissivity (-)

PropOpticalGround	=	struct('aveg',aveg_G,'abare',abare_G,'aimp',aimp_G,'albedo',albedo_G,...
						'eveg',eveg_G,'ebare',ebare_G,'eimp',eimp_G,'emissivity',emissivity_G);
					
% WALL OPTICAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
albedo_W		=	0.4;	% Wall surface albedo (-)
emissivity_W	=	0.95;	% Wall emissivity (-)

PropOpticalWall	=	struct('albedo',albedo_W,'emissivity',emissivity_W);

% TREE OPTICAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
albedo_T		=	0.2;	% Tree albedo (-)
emissivity_T	=	1 - exp(-(LAI_T+SAI_T));	% Tree emissivity (-)  

PropOpticalTree	=	struct('albedo',albedo_T,'emissivity',emissivity_T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THERMAL PROPERTIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROOF THERMAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lan_dry_imp_R	=	0.67;		% Thermal conductivity dry solid [W/m K]
cv_s_imp_R		=	1.0*10^6;	% Volumetric heat capacity solid [J/m^3 K]

ParThermalRoof	=	struct('lan_dry_imp',lan_dry_imp_R,'cv_s_imp',cv_s_imp_R);

% GROUND THERMAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lan_dry_imp_G	=	1.2;		% Thermal conductivity dry solid [W/m K]
cv_s_imp_G		=	1.5*10^6;	% Volumetric heat capacity solid [J/m^3 K]

ParThermalGround	=	struct('lan_dry_imp',lan_dry_imp_G,'cv_s_imp',cv_s_imp_G);

% WALL THERMAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lan_dry_W	=	0.67;			% Thermal conductivity dry solid [W/m K]
cv_s_W		=	1.0*10^6;		% Volumetric heat capacity solid [J/m^3 K]

ParThermalWall	=	struct('lan_dry',lan_dry_W,'cv_s',cv_s_W);

% VEGETATION THERMAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is not used in the model yet
Cthermal_leaf	=	640;	% [J m-2 K-1] Heat capacity per single leaf area based on Kitaya et al. 2003, Ryu et al. 2016

ParThermalTree	=	struct('Cthermal_leaf',Cthermal_leaf);            
					
					
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOLID LAYER DISCRETIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Layer thickness
dz1_R	=	0.1;       % Thickness of first roof layer [m]
dz2_R	=	0.1;       % Thickness of second roof layer [m]

% Soil layer discretization 
Zs_R	=	[ 0 10 20 50 100];	% soil layer discretization [mm]
ms_R	=	length(Zs_R)-1;		% number of soil layers [-]

% Fix soil moisture to field capacity in certain layers -> instead of
% irrigation at the top
FixSM_R             =   1; % 1=yes, 0, no
FixSM_LayerStart_R  =   1; % First layer with fixed soil moisture
FixSM_LayerEnd_R    =   ms_R; % Last layer with fixed soil moisture

% GROUND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soil layer discretization 
Zs_G	=	[0 10 20 50 100 150 200 300 400 600 800 1000 1500 2000];	% soil layer discretization [mm]
ms_G	=	length(Zs_G)-1;		% number of soil layers [-]

% Fix soil moisture to field capacity in certain layers -> instead of
% irrigation at the top
FixSM_G             =   1; % 1=yes, 0, no
FixSM_LayerStart_G  =   11; % First layer with fixed soil moisture
FixSM_LayerEnd_G    =   ms_G; % Last layer with fixed soil moisture

% WALL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Layer thickness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dz1_W	=	0.1;       % Thickness of first wall layer [m]
dz2_W	=	0.1;       % Thickness of second wall layer [m]

WallLayers	=	struct('dz1_wall',dz1_W,'dz2_wall',dz2_W);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTERCEPTION AND SOIL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soil classification: Sandy Loam
Pcla_R	=	0.20;    % Fraction of clay in the soil [-]
Psan_R	=	0.40;    % Fraction of sand in the soil [-]
Porg_R	=	0.025;   % Fraction of organic material in the soil [-]

% Interception and soil parameters
In_max_imp_R	=	0.25; 	% Maxiumum interception capacity of roof impervious area [mm]
In_max_ground_R	=	10;		% Maxiumum interception capacity of ground under roof vegetation [mm]
Sp_In_R			=	0.2;	% specific water retained by a vegetated surface [mm m^2 VEG area m^-2 plant area]

Kimp_R	=	0;				% Hydraulic conductivity of impervious area [mm/h]
Kfc_R	=	0.2;			% Conductivity at field capacity [mm/h]
Phy_R	=	10000;			% Suction at the residual/hygroscopic water content [kPa]
Kbot_R	=	NaN;			% [mm/h] Conductivity at the bedrock layer
SPAR_R	=	2;				% SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
% Do not use 1-VanGenuchten at the moment as very high soil water potential
% when dry

ParSoilRoof	=	struct('Zs',Zs_R,'ms',ms_R,'dz1',dz1_R,'dz2',dz2_R,...
                'FixSM_R',FixSM_R,'FixSM_LayerStart_R',FixSM_LayerStart_R,'FixSM_LayerEnd_R',FixSM_LayerEnd_R,...
				'In_max_imp',In_max_imp_R,'In_max_ground',In_max_ground_R,...
				'Sp_In',Sp_In_R,'Kimp',Kimp_R,'Kfc',Kfc_R,'Phy',Phy_R,'SPAR',SPAR_R,...
				'Kbot',Kbot_R,'Pcla',Pcla_R,'Psan',Psan_R,'Porg',Porg_R);

% GROUND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soil classification: Sandy Loam
Pcla_G	=	0.20;	% Fraction of clay in the soil [-]
Psan_G	=	0.40;	% Fraction of sand in the soil [-]
Porg_G	=	0.025;	% Fraction of organic material in the soil [-]

% Interception and soil parameters
In_max_imp_G		=	0.5; % 0.2;	% Maxiumum interception capacity of impervious ground area [mm]
In_max_underveg_G	=	10;		% Maxiumum interception capacity of vegetated ground area [mm]
In_max_bare_G		=	10;		% Maxiumum interception capacity of bare ground area [mm]
Sp_In_G				=	0.2;	% specific water retained by a vegetated surface on the ground [mm m^2 VEG area m^-2 plant area]

Kimp_G	=	0.001;		% Hydraulic conductivity of impervious area [mm/h]

Kfc_G	=	0.2;		% Conductivity at field capacity [mm/h]
Phy_G	=	10000;		% Suction at the residual/hygroscopic water content [kPa]
SPAR_G	=	2;			% SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
% Do not use 1-VanGenuchten at the moment as very high soil water potential
% when dry
Kbot_G	=	NaN;		% [mm/h] Conductivity at the bedrock layer

ParSoilGround	=	struct('Zs',Zs_G,'ms',ms_G,'In_max_imp',In_max_imp_G,...
                    'FixSM_G',FixSM_G,'FixSM_LayerStart_G',FixSM_LayerStart_G,'FixSM_LayerEnd_G',FixSM_LayerEnd_G,...
					'In_max_underveg',In_max_underveg_G,'In_max_bare',In_max_bare_G,...
					'Sp_In',Sp_In_G,'Kimp',Kimp_G,'Kfc',Kfc_G,'Phy',Phy_G,'SPAR',SPAR_G,...
					'Kbot',Kbot_G,'Pcla',Pcla_G,'Psan',Psan_G,'Porg',Porg_G);

% TREES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interception tree
Sp_In_T		=	0.2;	% specific water retained by the tree [mm m^2 VEG area m^-2 plant area]

ParInterceptionTree	=	struct('Sp_In',Sp_In_T);


% PERSON for MRT calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PositionPx		=	Gemeotry_m.Width_canyon/2;	% [m] position within canyon
PositionPz		=	1.1;						% [m] height of centre of person, usually choose 1.1 m

% PersonWidth and PersonHeight are not used at the moment
PersonWidth		=	0.06/2;		% [-] horizontal radius of ellipse describing person (=hip width / 2)
PersonHeight	=	0.22/2;		% [-] Vertical radius of ellipse describing person (= height / 2)

% Automatic wind speed calculation at user speficied height
HeightWind		=	1.1;		% [m] height for wind speed to calculate OTC

Person		=	struct('PositionPx',PositionPx,'PositionPz',PositionPz,...
	'PersonWidth',PersonWidth,'PersonHeight',PersonHeight,'HeightWind',HeightWind);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check feasibility of input parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_functions.GeometryCheck(Gemeotry_m,ParTree,Person,MeteoData.Zatm);
data_functions.InputParameterCheck(FractionsRoof,FractionsGround,ParVegRoof,ParVegGround,ParVegTree);


