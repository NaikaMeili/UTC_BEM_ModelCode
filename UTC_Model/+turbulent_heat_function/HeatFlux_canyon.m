function[Hcanyon,LEcanyon,ra_canyon,ra_orig,fconv,HumidityCan]=HeatFlux_canyon(TemperatureC,Gemeotry_m,MeteoData,ParVegTree,ParTree,fconvPreCalc,fconv)

% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hcanyon			=	sensible heat of canyon air (W/m^2)
% Ecanyon			=	latent heat of canyon air (x)

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zatm				=	Atmospheric reference height [m]
% Tatm				=	Air Temperature at atmospheric reference level [K]
% Uatm				=	Wind speed at atmospheric reference level [m/s]
% Pre				=	air pressure [Pa]. Carefull, used to be [hPa - mbar] in T&C
% T_canyon			=	Temperature of canyon air [K]
% Z0_canyon			=	roughness length of canyon [m]
% h_disp_canyon		=	displacement height of canyon [m]
% cp_atm			=	specific heat air  [J/kg K]
% rho_atm			=	dry air density at atmosphere [kg/m^3]
% q_canyon			=	Specific humidity of canyon []
% q_atm				=	Specifc humidity of air at reference height []

T_canyon		=	TemperatureC(1,9);
q_canyon		=	TemperatureC(1,10);
H				=	Gemeotry_m.Height_canyon;
W				=	Gemeotry_m.Width_canyon;
Wroof			=	Gemeotry_m.Width_roof;
Htree			=	Gemeotry_m.Height_tree;
R_tree			=	Gemeotry_m.Radius_tree;
Hcan_max		=	Gemeotry_m.Hcan_max;
Hcan_std		=	Gemeotry_m.Hcan_std;
Kopt			=	ParVegTree.Kopt;
LAI_t			=	ParVegTree.LAI;
trees			=	ParTree.trees;
Zatm			=	MeteoData.Zatm;
Tatm			=	MeteoData.Tatm;
Uatm			=	MeteoData.Uatm;
Pre				=	MeteoData.Pre;
q_atm			=	MeteoData.q_atm;
ea				=	MeteoData.ea;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensible heat flux canyon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cp_atm		=	1005+(((Tatm-273.15)+23.15)^2)/3364;		% specific heat air  [J/kg K]
rho_atm		=	(Pre/(287.04*Tatm))*(1-(ea/Pre)*(1-0.622));	% dry air density at atmosphere [kg/m^3]
% q_atm		=	0.622*ea/(Pre-0.378*ea);					% Specifc humidity of air at reference height []
L_heat		=	1000*(2501.3 - 2.361*(Tatm-273.15));		% Latent heat vaporization/condensaition [J/kg]

% % Vapor pressure saturation and specific humidity at esat
esat_T_canyon	=	611*exp(17.27*(T_canyon-273.16)/(237.3+(T_canyon-273.16)));	% vapor pressure saturation at T_canyon [Pa]
qsat_T_canyon	=	(0.622*esat_T_canyon)/(Pre-0.378*esat_T_canyon);			% Saturated specific humidity at T_canyon []
e_T_canyon		=	q_canyon*Pre/(0.622+0.378*q_canyon);
rel_hum_canyon	=	e_T_canyon/esat_T_canyon;

% Calculation of structural parameters and wind profile in the city
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dcan,zomcan,~,~,~,~]=resistance_functions.WindProfile_Canyon(...
	H,Htree,R_tree,W,Wroof,Kopt,LAI_t,Zatm,Uatm,H,trees,1.5,0.01,Hcan_max,Hcan_std);

zom_town	=	zomcan;			% Momentum roughness length of canyon, calculated according to McDonald 1998
zoh_town	=	zom_town/10;	% Heat roughness length of canyon

% Aerodynamic resistance calcuculated with Monin Obukhov similarity theory
[ra]=resistance_functions.Aerodynamic_Resistence(Tatm-273.15,T_canyon-273.15,Pre/100,Zatm,dcan,zom_town,zoh_town,Uatm,ea,e_T_canyon);

ra_orig = ra;

% Include enhancement term for aerodynamic resistance according to Pleim et
% al. 2007 accounting for non-local transport through large eddies
if fconvPreCalc ==1
    if T_canyon-Tatm > 0.1
        ra_enhanced = ra.*(1-fconv);
    else
       ra_enhanced = ra; 
    end
else
    if T_canyon-Tatm > 0.1
        hPBL = 1000;
        [fconv,ra_enhanced,~, LAN]=resistance_functions.EnhancementFactorRaPleim(ra,zom_town,zoh_town,dcan,Zatm,Uatm,hPBL);
        ra_enhanced = ra.*(1-fconv);
    else
       ra_enhanced = ra; 
 
    end
end

ra = ra_enhanced;

% Calculate heat fluxes
ra_canyon	=	ra;
Hcanyon    =	cp_atm*rho_atm*(T_canyon-Tatm)./ra;
Ecanyon	    =	rho_atm*(q_canyon-q_atm)./ra;
LEcanyon	=	L_heat*Ecanyon;


HumidityCan.CanyonRelative		=	rel_hum_canyon;
HumidityCan.CanyonSpecific		=	q_canyon;
HumidityCan.CanyonVapourPre		=	e_T_canyon;
HumidityCan.CanyonRelativeSat	=	1;
HumidityCan.CanyonSpecificSat	=	qsat_T_canyon;
HumidityCan.CanyonVapourPreSat	=	esat_T_canyon;

