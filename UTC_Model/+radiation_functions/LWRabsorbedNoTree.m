function[LWRin_nT,LWRout_nT,LWRabs_nT,LWREB_nT]...
         =LWRabsorbedNoTree(h_can,w_can,LWR,fgveg,fgbare,fgimp,ew,egveg,egbare,egimp,...
         Tgimp,Tgbare,Tgveg,Twsun,Twshade,ViewFactor)

% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LWRin_nT			=	total incoming longwave radiation to surface i [W/m^2]
% LWRout_nT			=	total outgoing longwave radiation of surface i [W/m^2]
% LWRabs_nT			=	total absorbed longwave radiation of surface i [W/m^2]
% LWREB_nT			=	Surface energy balance of each surface i [W/m^2]

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h_can			=	building height [-]
% w_can			=	ground width [-]
% fgveg			=	Partitioning ground vegetation [-]
% fgbare		=	Partitioning ground bare [-]
% fgimp			=	Partitioning ground impervious [-]
% egveg			=	ground emissivity (-)
% egbare		=	ground bare emissivity (-)
% egimp			=	ground impervious emissivity (-)
% erveg			=	Roof vegetation emissivity (-)
% erimp			=	Roof impervious emissivity (-)
% LWR			=	Atmospheric longwave radiation W per m^2 of horizontal surface [W/m^2]
% Tgimp			=	Temperature of the ground impervious [K]
% Tgbare		=	Temperature of the ground bare [K]
% Tgveg			=	Temperature of the ground vegetated [K]
% Twsun			=	Temperature of the wall sun [K]
% Twshade		=	Temperature of the wall shade [K]

F_gs_nT		=	ViewFactor.F_gs_nT;
F_gw_nT		=	ViewFactor.F_gw_nT;
F_ww_nT		=	ViewFactor.F_ww_nT;
F_wg_nT		=	ViewFactor.F_wg_nT;
F_ws_nT		=	ViewFactor.F_ws_nT;
F_sg_nT		=	ViewFactor.F_sg_nT;
F_sw_nT		=	ViewFactor.F_sw_nT;

% normalized surface areas
A_s		=	w_can;
A_g		=	w_can;
A_w		=	h_can;
bolzm	=	5.67*10^(-8);	% Stefan-Boltzmann constant [W*m^-2*K-4], Lecture notes hydrology 2

% Check if view factors add up to 1
% [F_gs_nT,F_gt_nT,F_gw_nT,F_ww_nT,F_wt_nT,F_wg_nT,F_ws_nT,F_ts_nT,F_tw_nT,F_tt_nT,F_tg_nT,F_sg_nT,F_sw_nT,F_st_nT]=radiation_functions.ViewFactorsNoTreeAnalytical(h);

SVF(1)			=	F_gs_nT + 2*F_gw_nT;
SVF(2)			=	F_ww_nT + F_wg_nT + F_ws_nT;
SVF(3)			=	F_sg_nT + 2*F_sw_nT;

SVF2(1)			=	F_gs_nT + 2*F_ws_nT*h_can;
SVF2(2)			=	F_sg_nT + 2*F_wg_nT*h_can;
SVF2(3)			=	F_ww_nT + F_sw_nT/h_can + F_gw_nT/h_can;

for i=length(SVF)
	if SVF(i)<0.999 || SVF(i)>1.001
	disp('The view factors do not add up to 1 for a canyon with trees')
	end
end

% Calculations for infinite reflections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sequence in Vectors :
% Vegetated ground
% Bare ground
% Impervious ground
% Sunlit wall
% Shaded wall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for infinite reflections equation A*X=C
Cimp	=	fgimp>0;
Cbare	=	fgbare>0;
Cveg	=	fgveg>0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% View factor matrix to solve for infinite reflections equation
% Omega_i = Tij*Bj
Tij	=	[1,0,0, -(1-egveg)*F_gw_nT*Cveg, -(1-egveg)*F_gw_nT*Cveg, -(1-egveg)*F_gs_nT*Cveg;...
		0,1,0, -(1-egbare)*F_gw_nT*Cbare, -(1-egbare)*F_gw_nT*Cbare, -(1-egbare)*F_gs_nT*Cbare;...
		0,0,1, -(1-egimp)*F_gw_nT*Cimp, -(1-egimp)*F_gw_nT*Cimp, -(1-egimp)*F_gs_nT*Cimp;...
		-(1-ew)*F_wg_nT*fgveg*Cveg, -(1-ew)*F_wg_nT*fgbare*Cbare, -(1-ew)*F_wg_nT*fgimp*Cimp, 1, -(1-ew)*F_ww_nT, -(1-ew)*F_ws_nT;...
		-(1-ew)*F_wg_nT*fgveg*Cveg, -(1-ew)*F_wg_nT*fgbare*Cbare, -(1-ew)*F_wg_nT*fgimp*Cimp, -(1-ew)*F_ww_nT, 1, -(1-ew)*F_ws_nT;...
		0, 0, 0, 0, 0, 1];
	
% Emitted radiation per surface
Omega_i	=	[(egveg*bolzm*(Tgveg)^4*Cveg);...
			(egbare*bolzm*(Tgbare)^4*Cbare);...
			(egimp*bolzm*(Tgimp)^4*Cimp);...
			(ew*bolzm*(Twsun)^4);...
			(ew*bolzm*(Twshade)^4);...
			LWR];
		

% How to solve the set of equations
% The outgoing and emitted radiation should be the same
% B_i				=	[Bveg; Bbare; Bimp; Bwall; Bwall; Bsky];
% Omega_i			=	Tij*B_i
% Tij^-1*Omega_i	=	Tij^-1*Tij*B_i
% Tij^-1*Omega_i	=	B_i

% Outgoing radiation per surface
B_i		=	Tij^-1*Omega_i;	% Outgoing radiation [W/m^2] per m^2 surface area

if B_i(6,1)~=LWR
	disp('Incoming lonwave radiation and emitted longwave radiation from the sky after the matrix inversion are not equal')
end

% Incoming longwave radiation at each surface A_i
Tij2	=	[0, 0, 0, F_gw_nT*Cveg, F_gw_nT*Cveg, F_gs_nT*Cveg;...
			0, 0, 0, F_gw_nT*Cbare, F_gw_nT*Cbare, F_gs_nT*Cbare;...
			0, 0, 0, F_gw_nT*Cimp, F_gw_nT*Cimp, F_gs_nT*Cimp;...
			F_wg_nT*fgveg*Cveg, F_wg_nT*fgbare*Cbare, F_wg_nT*fgimp*Cimp, 0, F_ww_nT, F_ws_nT;...
			F_wg_nT*fgveg*Cveg, F_wg_nT*fgbare*Cbare, F_wg_nT*fgimp*Cimp, F_ww_nT, 0, F_ws_nT;...
			0, 0, 0, 0, 0, 0];

A_i		=	Tij2*B_i;
e_i		=	[egveg; egbare; egimp; ew; ew; 0];
A_i2	=	(B_i-Omega_i)./(1-e_i);
Qnet_i2	=	A_i-B_i;


% Absorbed longwave radiation (Harman et al 2004)
e_i				=	[egveg; egbare; egimp; ew; ew; 0];

Qnet_i			=	(e_i.*B_i - Omega_i)./(1-e_i);
Qnet_i(e_i==1)	=	A_i(e_i==1) - Omega_i(e_i==1);
Qnet_i(6)		=	0; % Assumption: The sky has a fixed emission of LWR. Hence, Qnet is 0.

% Difference between emitted and outgoing
% Delta_out		=	B_i-Omega_i;

% Assignment
LWRout_i		=	B_i;	% Outgoing radiation [W/m^2] per m^2 surface area
LWRemit_i		=	Omega_i;% Emitted radiation [W/m^2] per m^2 surface area
LWRin_i			=	A_i;	% Incoming radiation [W/m^2] per m^2 surface area
LWRnet_i		=	Qnet_i;	% Net absorbed radiation [W/m^2] per m^2 surface area

% Energy balance
LWRin_atm			=	LWR;

TotalLWRSurface_in	=	LWRin_i(1)*fgveg*A_g/A_g + LWRin_i(2)*fgbare*A_g/A_g + LWRin_i(3)*fgimp*A_g/A_g + ...
						LWRin_i(4)*A_w/A_g + LWRin_i(5)*A_w/A_g;

TotalLWRSurface_abs	=	LWRnet_i(1)*fgveg*A_g/A_g + LWRnet_i(2)*fgbare*A_g/A_g + LWRnet_i(3)*fgimp*A_g/A_g + ...
						LWRnet_i(4)*A_w/A_g + LWRnet_i(5)*A_w/A_g;
      
TotalLWRSurface_out	=	LWRout_i(1)*fgveg*A_g/A_s+LWRout_i(2)*fgbare*A_g/A_s+LWRout_i(3)*fgimp*A_g/A_s+...
						LWRout_i(4)*A_w/A_s+LWRout_i(5)*A_w/A_s;
	  
TotalLWRref_to_atm	=	LWRout_i(1)*F_sg_nT*fgveg + LWRout_i(2)*F_sg_nT*fgbare + LWRout_i(3)*F_sg_nT*fgimp + ...
						LWRout_i(4)*F_sw_nT + LWRout_i(5)*F_sw_nT;

EBSurface			=	TotalLWRSurface_in - TotalLWRSurface_abs - TotalLWRSurface_out;
EBCanyon			=	LWRin_atm - TotalLWRSurface_abs - TotalLWRref_to_atm;

% Energy balance
if abs(EBSurface)>=10^-6
	disp('EBSurface is not 0. Please check LWRabsorbedNoTrees.m')
end
if abs(EBCanyon)>=10^-6
	disp('EBCanyon is not 0. Please check LWRabsorbedNoTrees.m')
end

% Sequence in Vectors :
% Vegetated ground
% Bare ground
% Impervious ground
% Sunlit wall
% Shaded wall

% Longwave radiation by each surface per m^2 surface area
% Incoming longwave radiation
LWRin_nT						=	[];			
LWRin_nT.LWRinGroundImp			=	LWRin_i(3)*Cimp;
LWRin_nT.LWRinGroundBare		=	LWRin_i(2)*Cbare;
LWRin_nT.LWRinGroundVeg			=	LWRin_i(1)*Cveg;
LWRin_nT.LWRinTree				=	0;
LWRin_nT.LWRinWallSun			=	LWRin_i(4);
LWRin_nT.LWRinWallShade			=	LWRin_i(5);
LWRin_nT.LWRinTotalGround		=	fgveg*LWRin_i(1)+fgbare*LWRin_i(2)+fgimp*LWRin_i(3);
LWRin_nT.LWRinTotalCanyon		=	LWRin_i(1)*fgveg*A_g/A_g + LWRin_i(2)*fgbare*A_g/A_g + LWRin_i(3)*fgimp*A_g/A_g + ...
									LWRin_i(4)*A_w/A_g + LWRin_i(5)*A_w/A_g;
								
% Outgoing longwave radiation
LWRout_nT						=	[];			
LWRout_nT.LWRoutGroundImp		=	LWRout_i(3)*Cimp;
LWRout_nT.LWRoutGroundBare		=	LWRout_i(2)*Cbare;
LWRout_nT.LWRoutGroundVeg		=	LWRout_i(1)*Cveg;
LWRout_nT.LWRoutTree			=	0;
LWRout_nT.LWRoutWallSun			=	LWRout_i(4);
LWRout_nT.LWRoutWallShade		=	LWRout_i(5);
LWRout_nT.LWRoutTotalGround		=	fgveg*LWRout_i(1)+fgbare*LWRout_i(2)+fgimp*LWRout_i(3);
LWRout_nT.LWRoutTotalCanyon		=	LWRout_i(1)*fgveg*A_g/A_g + LWRout_i(2)*fgbare*A_g/A_g + LWRout_i(3)*fgimp*A_g/A_g + ...
									LWRout_i(4)*A_w/A_g + LWRout_i(5)*A_w/A_g;
								
% Absorbed longwave radiation
LWRabs_nT						=	[];			
LWRabs_nT.LWRabsGroundImp		=	LWRnet_i(3)*Cimp;
LWRabs_nT.LWRabsGroundBare		=	LWRnet_i(2)*Cbare;
LWRabs_nT.LWRabsGroundVeg		=	LWRnet_i(1)*Cveg;
LWRabs_nT.LWRabsTree			=	0;
LWRabs_nT.LWRabsWallSun			=	LWRnet_i(4);
LWRabs_nT.LWRabsWallShade		=	LWRnet_i(5);
LWRabs_nT.LWRabsTotalGround		=	fgveg*LWRnet_i(1)+fgbare*LWRnet_i(2)+fgimp*LWRnet_i(3);
LWRabs_nT.LWRabsTotalCanyon		=	LWRnet_i(1)*fgveg*A_g/A_g + LWRnet_i(2)*fgbare*A_g/A_g + LWRnet_i(3)*fgimp*A_g/A_g + ...
									LWRnet_i(4)*A_w/A_g + LWRnet_i(5)*A_w/A_g;
								
% Energy Balance of longwave radiation								
LWREB_nT						=	[];			
LWREB_nT.LWREBGroundImp			=	LWRin_nT.LWRinGroundImp - LWRout_nT.LWRoutGroundImp - LWRabs_nT.LWRabsGroundImp;
LWREB_nT.LWREBGroundBare		=	LWRin_nT.LWRinGroundBare - LWRout_nT.LWRoutGroundBare - LWRabs_nT.LWRabsGroundBare;
LWREB_nT.LWREBGroundVeg			=	LWRin_nT.LWRinGroundVeg - LWRout_nT.LWRoutGroundVeg - LWRabs_nT.LWRabsGroundVeg;
LWREB_nT.LWREBTree				=	0;
LWREB_nT.LWREBWallSun			=	LWRin_nT.LWRinWallSun-LWRout_nT.LWRoutWallSun - LWRabs_nT.LWRabsWallSun;
LWREB_nT.LWREBWallShade			=	LWRin_nT.LWRinWallShade-LWRout_nT.LWRoutWallShade - LWRabs_nT.LWRabsWallShade;
LWREB_nT.LWREBTotalGround		=	LWRin_nT.LWRinTotalGround-LWRout_nT.LWRoutTotalGround - LWRabs_nT.LWRabsTotalGround;
LWREB_nT.LWREBTotalCanyon		=	LWRin_nT.LWRinTotalCanyon-LWRout_nT.LWRoutTotalCanyon - LWRabs_nT.LWRabsTotalCanyon;

if abs(LWREB_nT.LWREBGroundImp)>=10^-6 
	disp('LWREB_nT.LWREBGroundImp is not 0. Please check LWRabsorbedNoTrees.m')
end
if abs(LWREB_nT.LWREBGroundBare)>=10^-6
	disp('LWREB_nT.LWREBGroundBare is not 0. Please check LWRabsorbedNoTrees.m')
end
if abs(LWREB_nT.LWREBGroundVeg	)>=10^-6
	disp('LWREB_nT.LWREBGroundVeg	 is not 0. Please check LWRabsorbedNoTrees.m')
end
if abs(LWREB_nT.LWREBWallSun)>=10^-6
	disp('LWREB_nT.LWREBWallSun is not 0. Please check LWRabsorbedNoTrees.m')
end
if abs(LWREB_nT.LWREBWallShade)>=10^-6
	disp('LWREB_nT.LWREBWallShade is not 0. Please check LWRabsorbedNoTrees.m')
end
if abs(LWREB_nT.LWREBTotalGround)>=10^-6
	disp('LWREB_nT.LWREBTotalGround is not 0. Please check LWRabsorbedNoTrees.m')
end
if abs(LWREB_nT.LWREBTotalCanyon)>=10^-6
	disp('LWREB_nT.LWREBTotalCanyon is not 0. Please check LWRabsorbedNoTrees.m')
end

