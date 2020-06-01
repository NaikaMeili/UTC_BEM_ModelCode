function[dVlat]=SOIL_MOISTURES_RICH_COMP_LAT3(t,Vlat,...
	dz,SPAR,Ks,Osat,Ohy,L,Pe,O33,alpVG,nVG,...
	Cimp,Cbare,Cveg,fimp,fbare,fveg,Wcan) 


% Definition of assignment
% Vlat, Ohy, Osat, 
% 1: impervious soil column
% 2: bare soil column
% 3: vegetated soil column

%%% Need only - Oi 
Olat	=	Vlat./dz + Ohy;
I1		=	Olat >= Osat -1e-5 ; 
I2		=	Olat <= Ohy + 1e-5 ; 
Olat(I1)=	Osat -1e-5;
Olat(I2)=	Ohy + 1e-5; 

% Hydraulic conductivity and soil water potential
[Ko,Po]=soil_functions.Conductivity_Suction(SPAR,Ks,Osat,Ohy,L,Pe,O33,alpVG,nVG,Olat);

% Lateral water re-distribution
a				=	15;		% Assumption: horizontal and vertical unsaturated conductivity is the same
dxsoil			=	1000;	% [mm] = 1 [m]

Kf_gimp			=	Ko(1);	%[mm/h]
Kf_gbare		=	Ko(2);	%[mm/h]
Kf_gveg			=	Ko(3);	%[mm/h]

Psi_soil_gimp	=	Po(1);	%[mm]
Psi_soil_gbare	=	Po(2);	%[mm]
Psi_soil_gveg	=	Po(3);	%[mm]

% Soil water potential [mm]. The higher the soil water potential the drier
% the soil. Hence, a minus to change flux direction.
Qlat_bare2imp	=	-a.*(Kf_gbare+Kf_gimp)./2.*(Psi_soil_gbare-Psi_soil_gimp)./dxsoil;
Qlat_veg2imp	=	-a.*(Kf_gveg+Kf_gimp)./2.*(Psi_soil_gveg-Psi_soil_gimp)./dxsoil;
Qlat_veg2bare	=	-a.*(Kf_gbare+Kf_gveg)./2.*(Psi_soil_gveg-Psi_soil_gbare)./dxsoil;
Qlat_imp2bare	=	-a.*(Kf_gbare+Kf_gimp)./2.*(Psi_soil_gimp-Psi_soil_gbare)./dxsoil;
Qlat_bare2veg	=	-a.*(Kf_gbare+Kf_gveg)./2.*(Psi_soil_gbare-Psi_soil_gveg)./dxsoil;
Qlat_imp2veg	=	-a.*(Kf_gveg+Kf_gimp)./2.*(Psi_soil_gimp-Psi_soil_gveg)./dxsoil;

% Transmissivity
Tveg2imp	=	Qlat_veg2imp.*dz;
Tbare2imp	=	Qlat_bare2imp.*dz;
Tveg2bare	=	Qlat_veg2bare.*dz;
Timp2bare	=	Qlat_imp2bare.*dz;
Tbare2veg	=	Qlat_bare2veg.*dz;
Timp2veg	=	Qlat_imp2veg.*dz;

T_totflux	=	sum(nansum([Tbare2veg; Tveg2bare; Tbare2imp;...
				Timp2bare; Tveg2imp; Timp2veg],1));
			
if T_totflux ~= 0
	disp('The lateral transmissivities do not add up to 0. Please check SOIL_MOISTURES_RICH_COMP_LAT.m')
end

% Flux lateral (incoming flux 2imp, 2bare, 2veg are positive, outgoing
% negative)
Qin_veg2imp		=	Tveg2imp./(fimp*1000*Wcan)*Cveg*Cimp;		
Qin_bare2imp	=	Tbare2imp./(fimp*1000*Wcan)*Cbare*Cimp;	

Qin_veg2bare	=	Tveg2bare./(fbare*1000*Wcan)*Cveg*Cbare;		
Qin_imp2bare	=	Timp2bare./(fbare*1000*Wcan)*Cimp*Cbare;	

Qin_bare2veg	=	Tbare2veg./(fveg*1000*Wcan)*Cbare*Cveg;
Qin_imp2veg		=	Timp2veg./(fveg*1000*Wcan)*Cimp*Cveg;

% Total flux incoming to one soil column from the two other soil columns  [mm/h]
Qin_imp	=	nansum([Qin_veg2imp; Qin_bare2imp],1);
Qin_bare=	nansum([Qin_veg2bare; Qin_imp2bare],1);
Qin_veg	=	nansum([Qin_bare2veg; Qin_imp2veg],1);

%%%%%%%%%%%%%%%%%% SOIL WATER BALANCE without other terms %%%%%%%%%%%%%%%%%
dVlat	=	[Qin_imp; Qin_bare; Qin_veg];

Test	=	dVlat - [Qin_imp; Qin_bare; Qin_veg];

if isnan(sum(dVlat))
	disp('NaN values in the dVall')
end

