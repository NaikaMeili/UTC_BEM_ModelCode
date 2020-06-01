function[Zs,dz,ms,Osat,Ohy,nVG,alpVG,Ks_Zs,L,Pe,O33,SPAR,EvL_Zs,Inf_Zs,...
	RfH_Zs,RfL_Zs,Zinf,Kbot,Slo_pot,Dz,aR,aTop,rsd,lan_dry,lan_s,cv_s]...
	=SoilParametersTotal(Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,Zs)

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rrootl	=	Root length index [m root / m^2 PFT]
% PsiL50	=	[MPa]  Water Potential at 50% loss conductivity
% PsiX50	=	Water potential at 50 of xylem hyd
% Pcla		=	Fraction of clay in the soil [-]
% Psan		=	Fraction of sand in the soil [-]
% Porg		=	Fraction of organic material in the soil [-]
% Kfc		=	Conductivity at field capacity [mm/h]
% Phy		=	Suction at the residual/hygroscopic water content [kPa]
% SPAR		=	SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
% Kbot		=	[mm/h] Conductivity at the bedrock layer
% CASE_ROOT	=	Type of Root Profile
% ZR95		=	Root depth 95 percentile [mm]
% ZR50		=	Root depth 50 percentile [mm]
% ZRmax		=	Maximum Root depth [mm]
% Zs		=	soil layer discretization [mm]

% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zs		=	soil layer discretization [mm]
% dz		=	Thickness of soil layers [mm]
% ms		=	Number of soil layers [-]
% Osat		=	Saturation moisture 0 kPa [-]
% Ohy		=	[] Hygroscopic Moisture Evaporation cessation 10000 kPa - 10 MPa - 1000 m  
% nVG		=	n parameter Van-Genuchten soil water retention curve [1/mm]
% alpVG		=	Alpha parameter Van-Genuchten soil water retention curve [1/mm]
% Ks_Zs		=	Hydraulic conductivity at saturation for each soil layer [mm/h]
% L			=	Slope of logaritimc tension-moisture curve 
% Pe		=	Tension at air antry (bubbling pressure) [kPa]
% O33		=	33 kPa Moisture 
% SPAR		=	SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
% EvL_Zs	=	[-] Fraction of evaporation depth in a specific soil layer
% Inf_Zs	=	[-] Fraction of infiltration depth in a specific soil layer
% Rf_Zs		=	[-] Root Fraction in a given soil layer
% Rrootl	=	Root length index [m root / m^2 PFT]
% PsiL50	=	[MPa]  Water Potential at 50% loss conductivity
% PsiX50	=	Water potential at 50 of xylem hyd
% Zinf		=	[mm] Depth of infiltration layer (=first layer)
% Kbot		=	[mm/h] Conductivity at the bedrock layer
% Slo_pot	=	[fraction dy/dx]
% Dz		=	Delta Depth Between First Middle Layer and soil surface [mm]
% aR		=	anisotropy ratio
% aTop		=	[mm] Ratio betweeen Area/ContourLenght
% rsd		=	Normal density dry soil [kg/m^3]
% lan_dry	=	Thermal conductivity dry soil [W / m K]
% lan_s		=	thermal conductivity solid soil [W / m K]
% cv_s		=	Volumetric heat capacity solid soil [J / m^3 K]

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ms			=	length(Zs)-1;
dz			=	diff(Zs);			% [mm]  Thickness of the Layers
Dz			=	zeros(1,ms);

for i = 1:ms
    if i>1
        Dz(i)	=	(dz(i)+ dz(i-1))/2;	% Delta Depth Between Middle Layer  [mm]
    else
        Dz(i)	=	dz(1)/2;				% Delta Depth Between First Middle Layer and soil surface [mm]
    end
end

Zdes		=	Zs(2)-Zs(1); % [mm] Depth of desorption, Depth of evaporation layer (=first layer)
Zinf		=	Zs(2)-Zs(1); % [mm] Depth of infiltration layer (=first layer)

[EvL_Zs]=soil_functions.Evaporation_layers(Zs,Zdes); % [-] Fraction of evaporation depth in a specific soil layer
[Inf_Zs]=soil_functions.Evaporation_layers(Zs,Zinf); % [-] Fraction of infiltration depth in a specific soil layer

Slo_pot		=	zeros(1,ms);	% [fraction dy/dx]
aR			=	1;				% anisotropy ratio
cellsize	=	1;
aTop		=	1000*cellsize^2./cellsize; % [mm] Ratio betweeen Area/ContourLenght

%% Calculation of total thermal capacity out of soil composition for force restore method %%
%%%%% Soil Parameters: characterization of soil parameters out of soil composition %%
[Osat,L,Pe,Ks,O33,rsd,lan_dry,lan_s,cv_s,~]=soil_functions.Soil_parameters(Psan,Pcla,Porg);

rsd			=	rsd*ones(1,ms);		% Normal density dry soil [kg/m^3]
lan_dry		=	lan_dry*ones(1,ms);	% Thermal conductivity dry soil [W / m K]
lan_s		=	lan_s*ones(1,ms);	% Thermal conductivity solid soil [W / m K]
cv_s		=	cv_s*ones(1,ms);	% Volumetric heat capacity solid soil [J / m^3 K]

% Alpha parameter Van-Genuchten soil water retention curve [1/mm] %
p			=	3+2/L;
m			=	2/(p-1); 
nVG			=	1/(1-m);
alpVG		=	(((-101.9368*Pe)*(2*p*(p-1))/(p+3))*((55.6+7.4*p+p^2)/(147.8+8.1*p+0.092*p^2)))^-1; %%[1/mm]%;

%%% Initializing variables for each soil layer %%
Osat		=	Osat*ones(1,ms);	% Water content at saturation, saturation moisture 0 kPa  [-]
L			=	L*ones(1,ms);		% Slope of logarithmic tension-moisture curve [-]
Pe			=	Pe*ones(1,ms);		% Tension at air antry (bubbling pressure) [kPa]
Ks_Zs		=	Ks*ones(1,ms);		% Hydraulic conductivity at saturation for each soil layer [mm/h]
O33			=	O33*ones(1,ms);		% Soil water content at -33 [kPa] of water potential
nVG			=	nVG*ones(1,ms);		% n parameter Van-Genuchten soil water retention curve [1/mm]
alpVG		=	alpVG*ones(1,ms);	% Alpha parameter Van-Genuchten soil water retention curve [1/mm]

% Ohy [] Hygroscopic Moisture Evaporation cessation 10000 kPa - 10 MPa - 1000 m
% Pss		=	NaN		Pss is used for the calculation of Oss which we do not use here
% PwP		=	NaN		PwP is used for the calculation of Owp which we do not use here
% If the input into the soil_parameterII function is lower/equal to 12,
% SPAR = 2 is used. SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
[~,~,~,Ohy]=soil_functions.Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,NaN,NaN,Phy);

% Root Distribution ZR_H -ZR_L ---  Root Fraction in a given soil layer [1...m]
if CASE_ROOT_H == CASE_ROOT_L
	CASE_ROOT = CASE_ROOT_H;
else
	CASE_ROOT = CASE_ROOT_H;
	disp('CASE_ROOT_H and CASE_ROOT_L are not the same. CASE_ROOT_H is taken for the calculation')
end

[RfH_Zs,RfL_Zs]=soil_functions.Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L);




