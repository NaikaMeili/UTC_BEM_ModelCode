function[SunPosition,MeteoData,HumidityAtm,Anthropogenic,location,ParCalculation]...
	=UEHMForcingData(MeteoDataRaw,itt,varargin)

% Input variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LWR_in			=	MeteoDataRaw.LWR_in;		%[W/m2]
SAB1_in			=	MeteoDataRaw.SAB1_in;		%[W/m2]
SAB2_in			=	MeteoDataRaw.SAB2_in;		%[W/m2]
SAD1_in			=	MeteoDataRaw.SAD1_in;		%[W/m2]
SAD2_in			=	MeteoDataRaw.SAD2_in;		%[W/m2]
T_atm			=	MeteoDataRaw.T_atm;			%[K]
windspeed_u		=	MeteoDataRaw.windspeed_u;	%[m/s]
pressure_atm	=	MeteoDataRaw.pressure_atm;	%[Pa]
rain			=	MeteoDataRaw.rain;			%[mm/h]
rel_humidity	=	MeteoDataRaw.rel_humidity;	%[-]
date_time		=	MeteoDataRaw.Date;

%% Location properties of the urban area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi				=	47.38;			% latitude positive north (degrees)
lambda			=	8.56;			% longitude positive east (degrees)
theta_canyon	=	deg2rad(180);	% canyon orientation (rad)
DeltaGMT		=	1;				% difference with Greenwich Meridian Time [h]
		
location		=	struct('phi',phi,'lambda',lambda,'theta_canyon',theta_canyon,...
					'DeltaGMT',DeltaGMT);

%% Calculate zenith angle
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Datam			=	datevec(date_time);
Datam			=	Datam(itt,:);

t_bef			=	0.5;
t_aft			=	0.5;

[h_S,~,zeta_S,~,~,~,~] = data_functions.SetSunVariables(Datam, DeltaGMT, lambda, phi, t_bef,t_aft);

theta_Z			=	pi/2-h_S;	% solar zenith angle

if theta_Z<= -pi/2 || theta_Z>= pi/2
    theta_Z		=	pi/2;
end

theta_n			=	zeta_S-theta_canyon;	% difference between solar azimuth angle and canyon orientation

SunPosition		=	struct('Datam',Datam,'t_bef',t_bef,'t_aft',t_aft,...
					'theta_Z',theta_Z,'theta_n',theta_n,'zeta_S',zeta_S);

%% Radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SW_dir		=	SAB1_in(itt,1) + SAB2_in(itt,1);	% Direct incoming shortwave radiation W per m^2 of horizontal surface
SW_diff		=	SAD1_in(itt,1) + SAD2_in(itt,1);	% Diffuse incoming shortwave radiation W per m^2 of horizontal surface
LWR			=	LWR_in(itt,1);						% Atmospheric longwave radiation W per m^2 of horizontal surface

if abs(cos(theta_Z))< 0.1
    SW_dir		=	0;
end

%% Meteorological data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Zatm			=	25;						% Atmospheric reference height [m]
Tatm			=	T_atm(itt,1);			% Air Temperature at atmospheric reference level [K]
Uatm			=	windspeed_u(itt,1);		% Wind speed at atmospheric reference level [m/s]
Uatm(Uatm==0)	=	0.01;					% WINDSPEED CANNOT BE 0 OTHERWISE THE LEAF BOUNDARY RESISTANCE FAILS
esat_Tatm		=	611*exp(17.27*(Tatm-273.16)/(237.3+(Tatm-273.16)));	% vapor pressure saturation at Tatm [Pa]
rel_hum			=	rel_humidity(itt,1);	% Relative humidity [-]
ea				=	esat_Tatm*rel_hum;		% vapor pressure [Pa]
Pre				=	pressure_atm(itt,1);	% air pressure [Pa]. Carefull, used to be [hPa - mbar] in T&C
q_atm			=	0.622*ea/(Pre-0.378*ea);% Specifc humidity of air at reference height []
qSat_atm		=	0.622*esat_Tatm/(Pre-0.378*esat_Tatm); % Saturation specific humidity
Catm_CO2		=	400;					% [ppm]-[umolCO2/mol] Atmospheric CO2 concentration 2017
Catm_O2			=	210000;					% [ppm] - [umolO2/mol] Intercellular Partial Pressure Oxygen
Rain			=	rain(itt,1);			% Precipiation [mm]
Time			=	date_time(itt,1);	

MeteoData		=	struct('SW_dir',SW_dir,'SW_diff',SW_diff,'LWR',LWR,'Zatm',Zatm,...
					'Tatm',Tatm,'Uatm',Uatm,'esat_Tatm',esat_Tatm,'rel_hum',rel_hum,...
					'ea',ea,'Pre',Pre,'q_atm',q_atm,'Catm_CO2',Catm_CO2,'Catm_O2',Catm_O2,...
					'Rain',Rain,'Time',Time);

HumidityAtm		=	struct('AtmRelative',rel_hum,'AtmSpecific',q_atm,'AtmVapourPre',ea,...
					'AtmRelativeSat',1,'AtmSpecificSat',qSat_atm,'AtmVapourPreSat',esat_Tatm);

				

%% ANTHROPOGENIC FACTORS
% Building intertior temperature & anthropogenic heat input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (Tatm-273.15)<18			% Minimum temperature when building heating is switched on
	Tb = 18+273.15;
elseif (Tatm-273.15)>30		% Maximum temperature when air condition is switched on
	Tb = 30+273.15;
else	
	Tb	=	Tatm;			% The rest of the time, interior building temperature is equal to the air temperature at forcing height
end

Qf_canyon	=	10;	% Anthropogenic heat input into the canyon air [W/m^2]
Qf_roof		=	0;	% Not included in the model yet: Anthropogenic heat input above the roof [W/m^2]

% Anthropogenic water
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This irrigates so no water stress of plants occur. Feel free to change
% how you like.
if itt == 1
	Waterf_canyonVeg	=	0;		% [mm/h] applied on the vegetated ground surface area
else
	if varargin{1,1}.SoilPotWGroundTot_H(itt-1) <= -0.1 || varargin{1,1}.SoilPotWGroundVeg_L(itt-1)<=-0.3
		Waterf_canyonVeg	=	1;	% [mm/h] applied on the vegetated ground surface area
	else
		Waterf_canyonVeg	=	0;	% [mm/h] applied on the vegetated ground surface area
	end
end

Waterf_canyonBare	=	0;% [mm/h] applied on the vegetated ground surface area
Waterf_roof			=	0;% [mm/h] applied on the vegetated ground surface area

Anthropogenic	=	struct('Tb',Tb,'Qf_canyon',Qf_canyon,'Qf_roof',Qf_roof,...
					'Waterf_canyonVeg',Waterf_canyonVeg,'Waterf_canyonBare',Waterf_canyonBare,'Waterf_roof',Waterf_roof);

%% GENERAL PARAMETERS FOR CALCULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time steps other than dth=1 and tds=3600 are not extensively tested in this version of the code.
% If issues occur with varying the time step, please contact me (Naika Meili)
dth		=	1;										% time step of calculation [h]
dts		=	3600;									% time step of calculation [s]
row		=	1000;									% density of water [kg/m^3]
cp_atm	=	1005+(((Tatm-273.15)+23.15)^2)/3364;	% specific heat air  [J/kg K]
rho_atm	=	(Pre/(287.04*Tatm))*(1-(ea/Pre)*(1-0.622));	% dry air density at atmosphere [kg/m^3]
					
ParCalculation	=	struct('dts',dts,'dth',dth,'row',row,'cp_atm',cp_atm,'rho_atm',rho_atm);

return


