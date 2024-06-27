function[dcan,zomcan,u_Hcan_max,u_Zp,w_Zp,alpha,RoughnessParameter]=...
	WindProfile_Canyon(Hcan,Htree,R_tree,Wcan,Wroof,Kopt,LAI_t,Zatm,uatm,Zp,trees,Zref_und,zom_und,...
	Hcan_max,Hcan_std)

% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dcan		=	Urban displacement height including trees [m]
% zomcan	=	Urban momentum roughness height including trees [m]
% u_Hcan	=	Wind speed at canyon height [m/s]
% u_Zpcan	=	Wind speed within canyon at height Zpcan [m/s]
% w_Zpcan	=	Vertical wind speed within canyon [m/s]

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hcan		=	canyon height [m]
% Htree		=	Tree height [m]
% R_tree	=	Tree radius [m]
% Wcan		=	Canyon width [m]
% Wroof		=	Roof width [m]
% Kopt		=	Optical transmission factor [-]
% LAI_t		=	Leaf area index of tree [-]	
% Zatm		=	Atmospheric reference height [m]
% uatm		=	Wind speed at atmospheric reference height [m/s]
% Zpcan		=	Height of interest within the canyon [m]
% trees		=	Presence of trees [0: No, 1: Yes]


if R_tree==0
    trees=0;
end

if trees==0
	Htree	=	0;
	R_tree	=	0;
	LAI_t	=	0;
end


% Displacement height and roughness length according to Kent et al 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters McDonald 1998
a		=	4.43;   % Best fit for staggered arrays , or a = 3.59 for square arrays
k		=	0.4;	% Von Karman costant
b		=	1.0;	% b=1, no incorporation for drag correction factors. Good fit for staggered arrays
CDb		=	1.2;	% nominal drag for cubical obstacles

% Plan area fraction of buildings and vegetation
Ap_build	=	Wroof;
Ap_tree		=	4*R_tree;
Ap_urb		=	Wcan+Wroof;

% Frontal area fraction of vegetation and buildings: assumtpion infinite
% urban canyon perpendicular to the wind direction (Length of building and
% plot equals infinity)
TreeCrownHeight = min(3, 0.5*Hcan);

Af_build_s	=	Hcan*Wroof;
Af_veg_s	=	4*R_tree*TreeCrownHeight;
Ap_urbArea  =	(Wcan+Wroof)*Wroof;

% Tree canopy transmittance (optical = P2D)
P2D			=	exp(-Kopt*LAI_t);
P3D			=	P2D^0.40;			% Guan et al. 2003
Pv			=	(-1.251*P3D^2+ + 0.489*P3D + 0.803)/CDb;	% Guan et al. 2000

% Plan area fraction of buildings and 
% Calculation of structural parameters and wind profile in the city
Lp_tot		=	(Ap_build + (1-P3D)*Ap_tree)/Ap_urb;
H_tot		=	(Hcan*Ap_build + (Htree+R_tree)*(1-P3D)*Ap_tree)/(Ap_build + (1-P3D)*Ap_tree);

% Urban discplacement height and roughness length with incorportation of trees (Kent 2017), (MacDonald
% 1998)
dcan_MacD	=	H_tot*(1+a^(-Lp_tot)*(Lp_tot-1));	% displacement height of canyon [m], eq. 23

Af_build=	Af_build_s;
Af_veg	=	Af_veg_s;

zomcan_MacD	=	H_tot*(1-dcan_MacD/H_tot)*exp(-(1/k^2*0.5*b*CDb*(1-dcan_MacD/H_tot)*(Af_build+Pv*Af_veg)/Ap_urbArea)^(-0.5));

% Inclusion of Kanda et al (2013): Height variability of buildings and
% vegetation int the formulaton of dcan and zomcan
if not(isnan(Hcan_max))&& not(isnan(Hcan_std))
	a0	=	1.29;
	b0	=	0.36;
	c0	=	-0.17;

	X			=	(Hcan_std+H_tot)/Hcan_max;
	dcan_Kand	=	(c0.*X.^2+(a0.*Lp_tot.^b0-c0).*X).*Hcan_max;

	a1	=	0.71;
	b1	=	20.21;
	c1	=	-0.77;

	Y			=	Lp_tot*Hcan_std/H_tot;
	zomcan_Kand	=	(b1*Y^2+c1*Y+a1)*zomcan_MacD;

	dcan	=	dcan_Kand;
	zomcan	=	zomcan_Kand;
	zohcan	=	zomcan_Kand/10;
	RoughnessParameter	=	{'Kand'};
	
	if (dcan_Kand+zomcan_Kand) >= (Hcan_max-0.5)
			dcan		=	dcan_MacD;
			zomcan		=	zomcan_MacD;
			zohcan		=	zomcan_MacD/10;
			Hcan_max	=	Hcan;
			RoughnessParameter	=	{'MacD'};
% 			disp('dcan+zomcan > Hcan_max according to Kanda et al. 2012. Hence, the parametrization of MacDonald et al. 1998 is used.')
	end
	
else
	dcan		=	dcan_MacD;
	zomcan		=	zomcan_MacD;
	zohcan		=	zomcan_MacD/10;
	Hcan_max	=	Hcan;
	RoughnessParameter	=	{'MacD'};
end

% us_atm	=	k.*uatm./log((Zatm-dcan)./zomcan);		% Friction Velocity Atmosphere [m/s]
% Re		=	zomcan*us_atm/(1.46*10^(-5));	% Kanda2006
% Re		=	10;	% Bruesart
% Kb1	=	2.46*Re^(1/4)-2;	% Bruesart parameter of 2.46
% Kb2	=	1.29*Re^(1/4)-2;	% Kanda2006 parameter of 1.29
% zohcan1	=	zomcan/exp(Kb1)
% zohcan2	=	zomcan/exp(Kb2)

% Calculation of wind profile above and in the canyon with a logarithmic
% and exponential wind profile.
us_atm		=	k.*uatm./log((Zatm-dcan)./zomcan);		% Friction Velocity Atmosphere [m/s]
u_Hcan_max	=	(us_atm./k).*log((Hcan_max-dcan)./zomcan);	% Wind Speed at canyon top (highest building) [m/s]
alpha		=	log(uatm./u_Hcan_max)./(Zatm./Hcan_max - 1);	% Attenuation Coefficient Canyon not corrected for presence of trees.
if Zp >= Hcan_max
	u_Zp	=	(us_atm./k).*log((Zp-dcan)./zomcan);
	w_Zp	=	0;	
elseif Zp <= Hcan_max && Zp >= Zref_und
	u_Zp	=	u_Hcan_max.*exp(-alpha.*(1-Zp./Hcan_max));
	w_Zp	=	0;	
elseif Zp <= Zref_und && Zp >= zom_und
	uref_und	=	u_Hcan_max.*exp(-alpha.*(1-Zref_und./Hcan_max));
	usref_und	=	k*uref_und/log(Zref_und/zom_und);
	u_Zp		=	(usref_und./k).*log(Zp./zom_und);
	w_Zp		=	0;	
else
	u_Zp	=	0;
	w_Zp	=	0;	
	disp('wind speed calculation height higher than reference height or lower than roughness length')
end



