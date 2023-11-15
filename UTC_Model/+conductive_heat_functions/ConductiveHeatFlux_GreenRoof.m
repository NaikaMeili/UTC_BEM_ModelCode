function[G1,G2,dS]=ConductiveHeatFlux_GreenRoof(TemperatureR,TemperatureB,TempVec_ittm,Anthropogenic,Owater,...
										ParVegRoof,ParSoilRoof,ParThermalRoof,ParCalculation,BEM_on)

%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Troof		=	Surface temperature of roof [K]
% Tint		=	Temperature in the roof [K]
% Tint_tm1	=	Temperature in the roof at the previous time step [K]
% Tb		=	Interior temperature in the building [K]
% Otm1		=	Soil moisture at previous time step [-]
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

%%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% G1		=	Heat flux from surface to concrete interior [W/m^2]
% G2		=	Heat flux from concrete interior to building interior [W/m^2]
% dS		=	Energy storage in the roof

Troof		=	TemperatureR(1,2);
Tint		=	TemperatureR(1,4);
Tint_tm1	=	TempVec_ittm.TRoofIntVeg;
if BEM_on ==1
    Tb      =	TemperatureB(1); % Ceiling temperature
else
    Tb      =   Anthropogenic.Tb;
end
Otm1		=	Owater.OwRoofSoilVeg;
Rrootl		=	ParVegRoof.Rrootl;
PsiL50		=	ParVegRoof.PsiL50;
PsiX50		=	ParVegRoof.PsiX50;
CASE_ROOT	=	ParVegRoof.CASE_ROOT;
ZR95		=	ParVegRoof.ZR95;
ZR50		=	ParVegRoof.ZR50;
ZRmax		=	ParVegRoof.ZRmax;
Pcla		=	ParSoilRoof.Pcla;
Psan		=	ParSoilRoof.Psan;
Porg		=	ParSoilRoof.Porg;
Kfc			=	ParSoilRoof.Kfc;
Phy			=	ParSoilRoof.Phy;
SPAR		=	ParSoilRoof.SPAR;
Kbot		=	ParSoilRoof.Kbot;
Zs			=	ParSoilRoof.Zs;
dz1			=	ParSoilRoof.dz1;
dz2			=	ParSoilRoof.dz2;
cv_s2		=	ParThermalRoof.cv_s_imp;
lan_dry2	=	ParThermalRoof.lan_dry_imp;
dts			=	ParCalculation.dts;


% Troof		=	Surface temperature of roof [K]
% Tint		=	Temperature in the roof [K]
% Tint_tm1	=	Temperature in the roof at the previous time step [K]
% Tb		=	Interior temperature in the building [K]


%% Calculation of soil parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,dz,ms,Osat,Ohy,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,rsd,lan_dry,lan_s,cv_s]...
	=soil_functions.SoilParametersTotal(Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT,CASE_ROOT,0,ZR95,0,ZR50,0,ZRmax,Zs);

Tdamptm1	=	Tb;
Tdp(1,:)	=	Tdamptm1*ones(1,ms);	% Soil temperature of the layer
Tdptm1		=	Tdp(1,:);				% Soil temperature of the layer at the previous time step

[lanS,cv_soil,~]=soil_functions.Soil_Thermal_properties(Tdptm1-273.15,rsd,lan_dry,lan_s,cv_s,Osat,Ohy,Otm1);

% Average soil parameters according to soil layer thickness
dz			=	dz./sum(dz,2);
lanS		=	dz*lanS';
cv_soil		=	dz*cv_soil';

% Average roof parameters according to roof thickness
dz_roof		=	[dz1, dz2];
cv_roof		=	[cv_soil, cv_s2];
dz_roof		=	dz_roof./sum(dz_roof,2);
cv_roof		=	dz_roof*cv_roof';

%% Computation of heat fluxes
G1			=	lanS*(Troof-Tint)/dz1;	% Soil Heat Flux [W/m^2];
G2			=	lan_dry2*(Tint-Tb)/dz2;
dS			=	cv_roof*(dz1+dz2)/dts*(Tint-Tint_tm1);


