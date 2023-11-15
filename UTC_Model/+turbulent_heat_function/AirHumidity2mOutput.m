function[DEi,Eimp_2m,Ebare_soil_2m,Eveg_int_2m,Eveg_soil_2m,TEveg_2m,Ecan_2m,...
	q2m,e_T2m,RH_T2m,qcan,e_Tcan,RH_Tcan...
	]=AirHumidity2mOutput(q2m,T2m,Timp,Tbare,Tveg,Tcan,qcan,...
	rap_can2m,rap_can2m_Inv,rb_L,alp_soil_bare,r_soil_bare,alp_soil_veg,r_soil_veg,rs_sun_L,rs_shd_L,...
	dw_L,Fsun_L,Fshd_L,FractionsGround,ParVegGround,...
	Eimp,Ebare,Eveg_int,Eveg_pond,Eveg_soil,TEveg,Pre,...
    Humidity_ittm,fconv,MeteoData,Gemeotry_m,rho_atm,Zp1,ParCalculation)


% Vapor pressure and specific humidity at saturation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
esat_Timp	=	611*exp(17.27*(Timp-273.16)/(237.3+(Timp-273.16)));	% vapor pressure saturation at Tground_imp [Pa]
qsat_Timp	=	(0.622*esat_Timp)/(Pre-0.378*esat_Timp);			% Saturated specific humidity at Tground_imp []

esat_Tbare	=	611*exp(17.27*(Tbare-273.16)/(237.3+(Tbare-273.16)));% vapor pressure saturation at Tground_bare [Pa]
qsat_Tbare	=	(0.622*esat_Tbare)/(Pre-0.378*esat_Tbare);			% Saturated specific humidity at Tground_bare []

esat_Tveg	=	611*exp(17.27*(Tveg-273.16)/(237.3+(Tveg-273.16)));	% vapor pressure saturation at Tground_veg [Pa]
qsat_Tveg	=	(0.622*esat_Tveg)/(Pre-0.378*esat_Tveg);			% Saturated specific humidity at Tground_veg []

esat_Tcan	=	611*exp(17.27*(Tcan-273.16)/(237.3+(Tcan-273.16)));	% vapor pressure saturation at T_canyon [Pa]
esat_T2m	=	611*exp(17.27*(T2m-273.16)/(237.3+(T2m-273.16)));	% vapor pressure saturation at T_canyon [Pa]

e_Tcan		=	qcan*Pre/(0.622+0.378*qcan);	% Vapour pressure canyon  [Pa]
RH_Tcan		=	e_Tcan/esat_Tcan;				% Relative humidity canyon [-]

e_T2m		=	q2m*Pre/(0.622+0.378*q2m);	% Vapour pressure canyon  [Pa]
RH_T2m		=	e_T2m/esat_T2m;				% Relative humidity canyon [-]


%--------------------------------------------------------------------------
if Tcan-MeteoData.Tatm > 0.1
    ra_enhanced = rap_can2m_Inv.*(1-fconv);
else
   ra_enhanced = rap_can2m_Inv; 
end

% Turbulent heat fluxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eimp_2m_pot			=	rho_atm*(qsat_Timp-q2m)./rap_can2m;
Eimp_2m				=	min(Eimp_2m_pot,Eimp);

Ebare_soil_2m_pot	=	rho_atm*(alp_soil_bare*qsat_Tbare-q2m)./(rap_can2m+r_soil_bare);
Ebare_soil_2m		=	min(Ebare_soil_2m_pot,Ebare);

Eveg_int_2m_pot		=	rho_atm*(qsat_Tveg-q2m)./(rb_L/((ParVegGround.LAI+ParVegGround.SAI)*dw_L)+rap_can2m);
Eveg_int_2m			=	min(Eveg_int,Eveg_int_2m_pot);
Eveg_soil_2m_pot	=	rho_atm*(alp_soil_veg*qsat_Tveg-q2m)./(rap_can2m+r_soil_veg);
Eveg_soil_2m		=	min(Eveg_pond+Eveg_soil,Eveg_soil_2m_pot);
TEveg_sun_2m_pot	=	rho_atm*(qsat_Tveg-q2m)./(rb_L/((ParVegGround.LAI)*Fsun_L*(1-dw_L))+rap_can2m+rs_sun_L/((ParVegGround.LAI)*Fsun_L*(1-dw_L)));
TEveg_shd_2m_pot	=	rho_atm*(qsat_Tveg-q2m)./(rb_L/((ParVegGround.LAI)*Fshd_L*(1-dw_L))+rap_can2m+rs_shd_L/((ParVegGround.LAI)*Fshd_L*(1-dw_L)));
TEveg_2m_pot		=	TEveg_sun_2m_pot + TEveg_shd_2m_pot;
TEveg_2m			=	min(TEveg_2m_pot,TEveg);

Ecan_2m				=	rho_atm*(q2m-qcan)./ra_enhanced;

Vcanyon = (Gemeotry_m.Width_canyon*min(2*Zp1/Gemeotry_m.Height_canyon,1)*Gemeotry_m.Height_canyon)/Gemeotry_m.Width_canyon;
dS_E_air = Vcanyon*rho_atm*(q2m-Humidity_ittm.q2m)/ParCalculation.dts; % (W/m), m^2*(kg/m^3)*(J/kg)*(kg/kg)/s

% Equation set up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eimp_2m				=	FractionsGround.fimp*Eimp_2m;
Ebare_2m			=	FractionsGround.fbare*Ebare_soil_2m;
Eveg_2m				=	FractionsGround.fveg*(Eveg_int_2m+Eveg_soil_2m+TEveg_2m);

DEi = Ecan_2m+dS_E_air-Eimp_2m-Ebare_2m-Eveg_2m;
