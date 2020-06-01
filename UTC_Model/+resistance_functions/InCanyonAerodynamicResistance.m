function[rap_can,rap_Zp1,rap_Zp1_In,rap_Zp2,rap_Zp2_In,rap_Zp3,rap_Zp3_In,u_Hcan,u_Zp1,u_Zp2,u_Zp3,uref_und,alpha]...
	=InCanyonAerodynamicResistance(uatm,Zatm,Ta,Ts,hcan_max,hcan,dcan,zomcan,Zref_und,zom_und,Zp1,Zp2,Zp3,Pre,ea,RoughnessParameter)
% [rap_can,rap_Zp1,rap_Zp1_In,rap_Zp2,rap_Zp2_In,rap_Zp3,rap_Zp3_In,u_Hcan,u_Zp1,u_Zp2,u_Zp3,uref_und]=...
% resistance_functions.InCanyonAerodynamicResistance(10,40,20,25,30,20,0.123*30,1.5,0.1,0:0.1:40,4,10);

% Input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uatm		=	Wind speed at atmospheric reference height [m/s]
% Zatm		=	Atmospheric reference height [m]
% Ta		=	Air temperature [degrees Celcius]
% Ts		=	Surface temperature [degrees Celcius]
% hcan		=	Urban canyon height [m]
% dcan		=	Urban canyon displacement height [m]
% zomcan	=	Urban roughness length [m]
% Zref_und	=	Reference height within the canyon where exponential wind
%				profile changes to logarithmic again [m], e.g. 1.5 m
% zom_und	=	Roughness length of ground surface within the canyon.
%				Displacement height is assumed to be zero.
% Zp1		=	Height within the canyon of interest (Zref_und<Zp1<hcan)
% Zp2		=	Height within the canyon of interest (Zref_und<Zp2<hcan)
% Zp3		=	Height within the canyon of interest (Zref_und<Zp3<dcan+zomcan)

% Output parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rap_can		=	Aerodynamic urban undercanopy resistance from zom_und to
%					the canyon displacement height plus canyon roughness length[s/m]
% rap_Zp1		=	Undercanopy resistance from the ground to height Zp1 in
%					the urban canyon [s/m]
% rap_Zp1_In	=	Inverse undercanopy resistance from height Zp1 to the
%					canyon displacement height [s/m]
% rap_Zp2		=	Undercanopy resistance from the ground to height Zp2 in
%					the urban canyon [s/m]
% rap_Zp2_In	=	Inverse undercanopy resistance from height Zp2 to the
%					canyon displacement height [s/m]
% rap_Zp3		=	Undercanopy resistance from the ground to height Zp3 in
%					the urban canyon [s/m]
% rap_Zp3_In	=	Inverse undercanopy resistance from height Zp3 to the
%					canyon displacement height [s/m]
% u_Hcan		=	Wind speed at canyon height [m/s]
% u_Zp1			=	Wind speed at height Zp1 [m/s]
% u_Zp2			=	Wind speed at Zp2 [m/s]
% u_Zp3			=	Wind speed at Zp3 [m/s]
% uref_und		=	Wind speed at reference height within the canyon [m/s] 

% Height check and if necessary adjustment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Zref_und(Zref_und>(dcan+zomcan))
	disp('The urban undercanopy reference height is too big or the canyon too shallow: Zref_und>(dcan+zomcan)')
	Zref_und(Zref_und>=(dcan+zomcan))	=	dcan+zomcan-0.01;
end

if zom_und(zom_und>Zref_und)
	disp('The ground surface roughness is bigger than the urban undercanopy reference height: zom_und>Zref_und')
	zom_und(zom_und>=Zref_und)	=	Zref_und-0.01;
end

if isequal(RoughnessParameter{1},'MacD')
	hcan_max	=	hcan;
end

Zu_prof					=	Zp1;

Zp1(Zp1>(dcan+zomcan))	=	(dcan+zomcan);
Zp1(Zp1<Zref_und)		=	Zref_und;
Zp2(Zp2>(dcan+zomcan))	=	(dcan+zomcan);
Zp2(Zp2<Zref_und)		=	Zref_und;
Zp3(Zp3>(dcan+zomcan))	=	(dcan+zomcan);
Zp3(Zp3<Zref_und)		=	Zref_und;


% Constants
g		=	9.81; %%[m/s2]
k		=	0.41; %% Von Karman Constant

% Wind speeds and undercanopy resistance at Hcan and Zref
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
us_atm		=	k.*uatm./log((Zatm-dcan)./zomcan);			% Friction Velocity Atmosphere [m/s]
u_Hcan		=	(us_atm./k).*log((hcan_max-dcan)./zomcan);		% Wind Speed at canyon top [m/s]
alpha		=	log(uatm./u_Hcan)./(Zatm./hcan_max -1);			% Attenuation Coefficient Canyon not corrected for presence of trees.
uref_und	=	u_Hcan.*exp(-alpha.*(1-Zref_und./hcan_max));	% Reference height within the canyon
Kh_can		=	k^2.*uatm.*(hcan_max-dcan)./log((Zatm-dcan)./zomcan); % Eddy diffusion coefficient at canyon height

% Undercanopy resistance of canyon
rap_can		=	hcan_max.*exp(alpha)./(Kh_can.*alpha).*(exp(-alpha.*(Zref_und./hcan_max))-exp(-alpha.*((dcan+zomcan)./hcan_max)))...
				+ 1./(k^2.*uref_und)*log(Zref_und./zom_und).^2;
		
%%% Stability correction
P_Ref	=	100000; %% [Pa] reference pressure 
Rd		=	287.05; %% [J/kg K] dry air gas constant 
cp		=	1005 + ((Ta +23.15)^2)/3364; %% specific heat air  [J/kg K]
Oa		=	(Ta +273.15)*(Pre/P_Ref)^-(Rd/cp); %% Potential temperature air [K]
Os		=	(Ts +273.15)*(Pre/P_Ref)^-(Rd/cp); %% Potential temperature surface [K]

% esatTs	=	611*exp(17.27*(Os-273.16)/(237.3+(Os-273.16)));	% vapor pressure saturation at T_canyon [Pa]
% Ova		=	Oa*Pre./(Pre - 0.378*ea); %% Virtual Potential temperature air [K]
% Ovs		=	Os*Pre./(Pre - 0.378*esatTs); %% Virtual Potential temperature
Ova		=	Oa;
Ovs		=	Os;

Ri2 = (g.*(Ova-Ovs).*Zref_und)./(uref_und.^2.*(0.5.*(Ova+Ovs))); %% [-]
Ri2(Ri2>0.16)=0.16; %% Max. Stability
if Ri2 < 0 %% unstable
	rap_can = rap_can./((1-5.*Ri2).^(3/4));
else %% Stable
	rap_can = rap_can./((1-5.*Ri2).^2);
end


% Wind speeds and undercanopy resistance at Zp1 within the canyon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_Zp1	=	u_Hcan.*exp(-alpha.*(1-Zp1./hcan_max)); %% Wind at canyon height Zp
rap_Zp1	=	hcan_max.*exp(alpha)./(Kh_can.*alpha).*(exp(-alpha.*(Zref_und./hcan_max))-exp(-alpha.*(Zp1./hcan_max)))...
			+ 1./(k^2.*uref_und).*log(Zref_und./zom_und).^2;

if Ri2 < 0 %% unstable
	rap_Zp1 = rap_Zp1./((1-5.*Ri2).^(3/4));
else %% Stable
	rap_Zp1 = rap_Zp1./((1-5.*Ri2).^2);
end

rap_Zp1_In	=	max(rap_can-rap_Zp1,0.1);
% rap_Zp1_In cannot be 0 as this lead to a division by 0. Hence, it is set
% to a small value (0.1)


% Wind speeds and undercanopy resistance at Zp2 within the canyon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_Zp2	=	u_Hcan.*exp(-alpha.*(1-Zp2./hcan_max)); %% Wind at canyon height Zp
rap_Zp2	=	hcan_max.*exp(alpha)./(Kh_can.*alpha).*(exp(-alpha.*(Zref_und./hcan_max))-exp(-alpha.*(Zp2./hcan_max)))...
			+ 1./(k^2.*uref_und).*log(Zref_und./zom_und).^2;

if Ri2 < 0 %% unstable
	rap_Zp2 = rap_Zp2./((1-5.*Ri2).^(3/4));
else %% Stable
	rap_Zp2 = rap_Zp2./((1-5.*Ri2).^2);
end

rap_Zp2_In	=	max(rap_can-rap_Zp2,0.1);

% Wind speeds and undercanopy resistance at Zp3 within the canyon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_Zp3	=	u_Hcan.*exp(-alpha.*(1-Zp3./hcan_max)); %% Wind at canyon height Zp
rap_Zp3	=	hcan_max.*exp(alpha)./(Kh_can.*alpha).*(exp(-alpha.*(Zref_und./hcan_max))-exp(-alpha.*(Zp3./hcan_max)))...
			+ 1./(k^2.*uref_und).*log(Zref_und./zom_und).^2;

if Ri2 < 0 %% unstable
	rap_Zp3 = rap_Zp3./((1-5.*Ri2).^(3/4));
else %% Stable
	rap_Zp3 = rap_Zp3./((1-5.*Ri2).^2);
end

rap_Zp3_In	=	max(rap_can-rap_Zp3,0.1);

% % Wind profile
% u_log1		=	(us_atm./k).*log((Zu_prof(Zu_prof>=hcan_max)-dcan)./zomcan);	% Wind Speed at canyon top [m/s]
% u_exp		=	u_Hcan.*exp(-alpha.*(1-Zu_prof(Zu_prof<=hcan_max&Zu_prof>=Zref_und)./hcan_max));		% Wind at canyon height Zp [m/s]
% usref_und	=	k*uref_und/log(Zref_und/zom_und);			% Friction Velocity Atmosphere [m/s]
% u_log2		=	(usref_und./k).*log(Zu_prof(Zu_prof<=Zref_und)./zom_und);		% Wind Speed at canyon ground [m/s]
% 
% figure
% plot(u_log1,Zu_prof(Zu_prof>=hcan_max),'DisplayName','Wind speed above canyon[m/s]')
% hold on
% plot(u_exp,Zu_prof(Zu_prof<=hcan_max&Zu_prof>=Zref_und),'DisplayName','Wind speed in canyon [m/s]')
% plot(u_log2,Zu_prof(Zu_prof<=Zref_und),'DisplayName','Wind speed at ground surface [m/s]')
% xlabel('Canyon height [m]')
% ylabel('mean wind speed u [m/s]')
% legend('show')
% 
% 
% figure
% plot(rap_Zp1,Zp1,'DisplayName','Undercanopy resistance from the ground to height Z [s/m]')
% % hold on
% % plot(rap_Zp1_In,Zp1,'DisplayName','Undercanopy resistance from canyon air to height z [s/m]')
% xlabel('Canyon height [m]')
% ylabel('Undercanopy resistance [s/m]')
% legend('show')


