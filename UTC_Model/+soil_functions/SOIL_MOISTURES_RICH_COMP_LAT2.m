function[dVlat]=SOIL_MOISTURES_RICH_COMP_LAT2(t,Vlat,...
	dz,SPAR,Ks,Osat,Ohy,L,Pe,O33,alpVG,nVG,C1,C2,f1,f2,Wcan) 

% Definition of assignment
% Vlat, Ohy, Osat, 

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

% Soil water potential [mm]. The higher the soil water potential the drier
% the soil. Hence, I put a minus to change flux direction.
Qlat_1to2	=	-a.*(Ko(1)+Ko(2))./2.*(Po(1)-Po(2))./dxsoil;
Qlat_2to1	=	-a.*(Ko(1)+Ko(2))./2.*(Po(2)-Po(1))./dxsoil;

% Transmissivity
T_1to2	=	Qlat_1to2.*dz;
T_2to1	=	Qlat_2to1.*dz;

T_totflux	=	sum(nansum([T_1to2; T_2to1],1));
			
if T_totflux ~= 0
	disp('The lateral transmissivities do not add up to 0. Please check SOIL_MOISTURES_RICH_COMP_LAT.m')
end

% Flux In
Qin_1to2	=	T_1to2./(f2*1000*Wcan)*C1*C2;
Qin_2to1	=	T_2to1./(f1*1000*Wcan)*C2*C1;	

% Total flux incoming to one soil column from the two other soil columns  [mm/h]
Qin_1		=	nansum(Qin_2to1,1);
Qin_2		=	nansum(Qin_1to2,1);

%%%%%%%%%%%%%%%%%% SOIL WATER BALANCE without other terms %%%%%%%%%%%%%%%%%
dVlat	=	[Qin_1; Qin_2];

Test	=	dVlat - [Qin_1; Qin_2];

if isnan(sum(dVlat))
	disp('NaN values in the dVall')
end

