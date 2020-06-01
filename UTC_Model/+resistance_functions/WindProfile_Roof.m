function[u_Zp,u_Hveg]=WindProfile_Roof(Hcan,Zatm,uatm,zom,disp_h,hveg,Zp)
% [u_Zp]=resistance_functions.WindProfile_Roof(10,20,5,0.01,0:0.1:10);

% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% u_Zp	=	Wind speed within at height Zp [m/s]

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hcan		=	canyon height [m]
% Zatm		=	Atmospheric reference height [m]
% uatm		=	Wind speed at atmospheric reference height [m/s]
% Zp		=	Height of interest [m]
% zom		=	Roughness length of roof cover [m]

% Calculation of wind profile above and in the canyon with a logarithmic
% wind profile.
k		=	0.4;	% Von Karman costant
Zatm	=	Zatm-Hcan;
us_atm	=	k.*uatm./log(Zatm./zom);	% Friction Velocity Atmosphere [m/s]
u_Zp	=	(us_atm./k).*log((Zp-disp_h)./zom);	% Wind Speed at height Zp [m/s]
u_Hveg	=	(us_atm./k).*log((hveg-disp_h)./zom);	% Wind Speed at vegetation top [m/s]



