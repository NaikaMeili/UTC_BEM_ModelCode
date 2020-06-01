function[LWRin_T,LWRout_T,LWRabs_T,LWREB_T]...
         =LWRabsorbedWithTrees(h_can,w_can,r_tree,LWR,fgveg,fgbare,fgimp,ew,et,egveg,egbare,egimp,...
		 Tgimp,Tgbare,Tgveg,Twsun,Twshade,Ttree,ViewFactor)
     
% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LWRin_T		=	total incoming longwave radiation to surface i [W/m^2]
% LWRout_T		=	total outgoing longwave radiation of surface i [W/m^2]
% LWRabs_T		=	total absorbed longwave radiation of surface i [W/m^2]
% LWREB_T		=	Surface energy balance of each surface i [W/m^2]

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h_can		=	building height [-]
% w_can		=	ground width [-]
% d_tree	=	location of trees in the canyon, tree-wall distance  [-]
% h_tree	=	height of trees, vertical level at the crown center  [-]
% r_tree	=	size of the tree crown, crown radius  [-]
% fgveg		=	Partitioning ground vegetation [-]
% fgbare	=	Partitioning ground bare [-]
% fgimp		=	Partitioning ground impervious [-]
% ew		=	Wall emissivity (-)
% et		=	tree emissivity (-)
% egveg		=	ground vegetation emissivity (-)
% egbare	=	ground bare emissivity (-)
% egimp		=	ground impervious emissivity (-)
% LWR		=	Atmospheric longwave radiation W per m^2 of horizontal surface
% Tgimp		=	Temperature of the ground impervious [K]
% Tgbare	=	Temperature of the ground bare [K]
% Tgveg		=	Temperature of the ground vegetated [K]
% Twsun		=	Temperature of the wall sun [K]
% Twshade	=	Temperature of the wall shade [K]
% Ttree		=	Temperature of the tree [K]

F_gs_T		=	ViewFactor.F_gs_T;
F_gt_T		=	ViewFactor.F_gt_T;
F_gw_T		=	ViewFactor.F_gw_T;
F_ww_T		=	ViewFactor.F_ww_T;
F_wt_T		=	ViewFactor.F_wt_T;
F_wg_T		=	ViewFactor.F_wg_T;
F_ws_T		=	ViewFactor.F_ws_T;
F_sg_T		=	ViewFactor.F_sg_T;
F_sw_T		=	ViewFactor.F_sw_T;
F_st_T		=	ViewFactor.F_st_T;
F_tg_T		=	ViewFactor.F_tg_T;
F_tw_T		=	ViewFactor.F_tw_T;
F_ts_T		=	ViewFactor.F_ts_T;
F_tt_T		=	ViewFactor.F_tt_T;

% % normalized surface areas
A_s		=	w_can;
A_g		=	w_can;
A_w		=	h_can;
A_t		=	2*2*pi*r_tree;        % There are 2 trees. Hence, the area of tree is twice a circle
bolzm	=	5.67*10^(-8);          % Stefan-Boltzmann constant [W*m^-2*K-4], Lecture notes hydrology 2

% Check if view factors add up to 1
SVF(1)			=	F_gs_T+F_gt_T+2*F_gw_T;
SVF(2)			=	F_ww_T+F_wt_T+F_wg_T+F_ws_T;
SVF(3)			=	F_sg_T+2*F_sw_T+F_st_T;
SVF(4)			=	F_ts_T+2*F_tw_T+F_tt_T+F_tg_T;

SVF2(1)			=	F_gs_T+2*F_ws_T+F_ts_T;
SVF2(2)			=	F_sg_T+2*F_wg_T+F_tg_T;
SVF2(3)			=	F_ww_T+F_sw_T+F_gw_T+F_tw_T;
SVF2(4)			=	F_gt_T+2*F_wt_T+F_tt_T+F_st_T;

for i=length(SVF)
	if SVF(i)<0.999 || SVF(i)>1.001
	disp('The view factors do not add up to 1 for a canyon with trees')
	end
end


% Calculations for infinite reflections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sequence in Vectors :
%
% Vegetated ground
% Bare ground
% Impervious ground
% Sunlit wall
% Shaded wall
% Trees
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for infinite reflections equation A*X=C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cimp	=	fgimp>0;
Cbare	=	fgbare>0;
Cveg	=	fgveg>0;

% View factor matrix to solve for infinite reflections equation
% Omega_i = Tij*Bj
Tij	=	[1,0,0, -(1-egveg)*F_gw_T*Cveg, -(1-egveg)*F_gw_T*Cveg, -(1-egveg)*F_gt_T*Cveg, -(1-egveg)*F_gs_T*Cveg;...
		0,1,0, -(1-egbare)*F_gw_T*Cbare, -(1-egbare)*F_gw_T*Cbare, -(1-egbare)*F_gt_T*Cbare, -(1-egbare)*F_gs_T*Cbare;...
		0,0,1, -(1-egimp)*F_gw_T*Cimp, -(1-egimp)*F_gw_T*Cimp, -(1-egimp)*F_gt_T*Cimp, -(1-egimp)*F_gs_T*Cimp;...
		-(1-ew)*F_wg_T*fgveg*Cveg, -(1-ew)*F_wg_T*fgbare*Cbare, -(1-ew)*F_wg_T*fgimp*Cimp, 1, -(1-ew)*F_ww_T, -(1-ew)*F_wt_T, -(1-ew)*F_ws_T;...
		-(1-ew)*F_wg_T*fgveg*Cveg, -(1-ew)*F_wg_T*fgbare*Cbare, -(1-ew)*F_wg_T*fgimp*Cimp, -(1-ew)*F_ww_T, 1, -(1-ew)*F_wt_T, -(1-ew)*F_ws_T;...
		-(1-et)*F_tg_T*fgveg*Cveg, -(1-et)*F_tg_T*fgbare*Cbare, -(1-et)*F_tg_T*fgimp*Cimp, -(1-et)*F_tw_T, -(1-et)*F_tw_T, 1-(1-et)*F_tt_T, -(1-et)*F_ts_T;...
		0, 0, 0, 0, 0, 0, 1];
	
% 1-((1-ew)*F_tt_T)
	
% Emitted radiation per surface
Omega_i	=	[(egveg*bolzm*(Tgveg)^4*Cveg);...
			(egbare*bolzm*(Tgbare)^4*Cbare);...
			(egimp*bolzm*(Tgimp)^4*Cimp);...
			(ew*bolzm*(Twsun)^4);...
			(ew*bolzm*(Twshade)^4);...
			(et*bolzm*(Ttree)^4);...
			LWR];
	
% How to solve the set of equations
% The outgoing and emitted radiation should be the same
% B_i				=	[Bveg; Bbare; Bimp; Bwall; Bwall; Bsky];
% Omega_i			=	Tij*B_i
% Tij^-1*Omega_i	=	Tij^-1*Tij*B_i
% Tij^-1*Omega_i	=	B_i

% Outgoing radiation per surface
B_i		=	Tij^-1*Omega_i;	% Outgoing radiation [W/m^2] per m^2 surface area

if B_i(7,1)~=LWR
	disp('Incoming lonwave radiation and emitted longwave radiation from the sky after the matrix inversion are not equal')
end

% Incoming longwave radiation at each surface A_i
Tij2=	[0, 0, 0, F_gw_T*Cveg, F_gw_T*Cveg, F_gt_T*Cveg, F_gs_T*Cveg;...
		0, 0, 0, F_gw_T*Cbare, F_gw_T*Cbare, F_gt_T*Cbare, F_gs_T*Cbare;...
		0, 0, 0, F_gw_T*Cimp, F_gw_T*Cimp, F_gt_T*Cimp, F_gs_T*Cimp;...
		F_wg_T*fgveg*Cveg, F_wg_T*fgbare*Cbare, F_wg_T*fgimp*Cimp, 0, F_ww_T, F_wt_T, F_ws_T;...
		F_wg_T*fgveg*Cveg, F_wg_T*fgbare*Cbare, F_wg_T*fgimp*Cimp, F_ww_T, 0, F_wt_T, F_ws_T;...
		F_tg_T*fgveg*Cveg, F_tg_T*fgbare*Cbare, F_tg_T*fgimp*Cimp, F_tw_T, F_tw_T, F_tt_T, F_ts_T;...
		0, 0, 0, 0, 0, 0, 0];

%F_tt_T
	
A_i		=	Tij2*B_i;

e_i		=	[egveg; egbare; egimp; ew; ew; et; 0];
A_i2	=	(B_i-Omega_i)./(1-e_i);
Qnet_i2	=	A_i-B_i;


% Absorbed longwave radiation (Harman et al 2004)
e_i				=	[egveg; egbare; egimp; ew; ew; et; 1];

Qnet_i			=	(e_i.*B_i - Omega_i)./(1-e_i);
Qnet_i(e_i==1)	=	A_i(e_i==1) - Omega_i(e_i==1);
Qnet_i(7)		=	0; % Assumption: The sky has a fixed emission of LWR. Hence, Qnet is 0.

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
						LWRin_i(4)*A_w/A_g + LWRin_i(5)*A_w/A_g + LWRin_i(6)*A_t/A_g;

TotalLWRSurface_abs	=	LWRnet_i(1)*fgveg*A_g/A_g + LWRnet_i(2)*fgbare*A_g/A_g + LWRnet_i(3)*fgimp*A_g/A_g + ...
						LWRnet_i(4)*A_w/A_g + LWRnet_i(5)*A_w/A_g + LWRnet_i(6)*A_t/A_g;
      
TotalLWRSurface_out	=	LWRout_i(1)*fgveg*A_g/A_s+LWRout_i(2)*fgbare*A_g/A_s+LWRout_i(3)*fgimp*A_g/A_s+...
						LWRout_i(4)*A_w/A_s+LWRout_i(5)*A_w/A_s+LWRout_i(6)*A_t/A_s;
	  
TotalLWRref_to_atm	=	LWRout_i(1)*F_sg_T*fgveg + LWRout_i(2)*F_sg_T*fgbare + LWRout_i(3)*F_sg_T*fgimp + ...
						LWRout_i(4)*F_sw_T + LWRout_i(5)*F_sw_T + LWRout_i(6)*F_st_T;

EBSurface			=	TotalLWRSurface_in - TotalLWRSurface_abs - TotalLWRSurface_out;
EBCanyon			=	LWRin_atm - TotalLWRSurface_abs - TotalLWRref_to_atm;

% Energy balance
if abs(EBSurface)>=10^-6
	disp('EBSurface LWR is not 0. Please check LWRabsorbedWithTrees.m')
end
if abs(EBCanyon)>=10^-6
	disp('EBCanyon LWR is not 0. Please check LWRabsorbedWithTrees.m')
end

if abs(EBCanyon)>=10^-6
	msgbox('EBCanyon LWR is not 0, calculate the appropriate view factors', 'Error','error');
	error('EBCanyon LWR is not 0, calculate the appropriate view factors')	
end	

% Sequence in Vectors :
% Vegetated ground
% Bare ground
% Impervious ground
% Sunlit wall
% Shaded wall
% Trees
% Longwave radiation by each surface per m^2 surface area
% Incoming longwave radiation
LWRin_T							=	[];			
LWRin_T.LWRinGroundImp			=	LWRin_i(3)*Cimp;
LWRin_T.LWRinGroundBare			=	LWRin_i(2)*Cbare;
LWRin_T.LWRinGroundVeg			=	LWRin_i(1)*Cveg;
LWRin_T.LWRinTree				=	LWRin_i(6);
LWRin_T.LWRinWallSun			=	LWRin_i(4);
LWRin_T.LWRinWallShade			=	LWRin_i(5);
LWRin_T.LWRinTotalGround		=	fgveg*LWRin_i(1)+fgbare*LWRin_i(2)+fgimp*LWRin_i(3);
LWRin_T.LWRinTotalCanyon		=	LWRin_i(1)*fgveg*A_g/A_g + LWRin_i(2)*fgbare*A_g/A_g + LWRin_i(3)*fgimp*A_g/A_g + ...
									LWRin_i(4)*A_w/A_g + LWRin_i(5)*A_w/A_g + LWRin_i(6)*A_t/A_g;						
% Outgoing longwave radiation
LWRout_T						=	[];			
LWRout_T.LWRoutGroundImp		=	LWRout_i(3)*Cimp;
LWRout_T.LWRoutGroundBare		=	LWRout_i(2)*Cbare;
LWRout_T.LWRoutGroundVeg		=	LWRout_i(1)*Cveg;
LWRout_T.LWRoutTree				=	LWRout_i(6);
LWRout_T.LWRoutWallSun			=	LWRout_i(4);
LWRout_T.LWRoutWallShade		=	LWRout_i(5);
LWRout_T.LWRoutTotalGround		=	fgveg*LWRout_i(1)+fgbare*LWRout_i(2)+fgimp*LWRout_i(3);
LWRout_T.LWRoutTotalCanyon		=	LWRout_i(1)*fgveg*A_g/A_g + LWRout_i(2)*fgbare*A_g/A_g + LWRout_i(3)*fgimp*A_g/A_g + ...
									LWRout_i(4)*A_w/A_g + LWRout_i(5)*A_w/A_g + LWRout_i(6)*A_t/A_g;						
% Absorbed longwave radiation
LWRabs_T						=	[];			
LWRabs_T.LWRabsGroundImp		=	LWRnet_i(3)*Cimp;
LWRabs_T.LWRabsGroundBare		=	LWRnet_i(2)*Cbare;
LWRabs_T.LWRabsGroundVeg		=	LWRnet_i(1)*Cveg;
LWRabs_T.LWRabsTree				=	LWRnet_i(6);
LWRabs_T.LWRabsWallSun			=	LWRnet_i(4);
LWRabs_T.LWRabsWallShade		=	LWRnet_i(5);
LWRabs_T.LWRabsTotalGround		=	fgveg*LWRnet_i(1)+fgbare*LWRnet_i(2)+fgimp*LWRnet_i(3);
LWRabs_T.LWRabsTotalCanyon		=	LWRnet_i(1)*fgveg*A_g/A_g + LWRnet_i(2)*fgbare*A_g/A_g + LWRnet_i(3)*fgimp*A_g/A_g + ...
									LWRnet_i(4)*A_w/A_g + LWRnet_i(5)*A_w/A_g + LWRnet_i(6)*A_t/A_g;
% Energy Balance of longwave radiation								
LWREB_T							=	[];			
LWREB_T.LWREBGroundImp			=	LWRin_T.LWRinGroundImp - LWRout_T.LWRoutGroundImp - LWRabs_T.LWRabsGroundImp;
LWREB_T.LWREBGroundBare			=	LWRin_T.LWRinGroundBare - LWRout_T.LWRoutGroundBare - LWRabs_T.LWRabsGroundBare;
LWREB_T.LWREBGroundVeg			=	LWRin_T.LWRinGroundVeg - LWRout_T.LWRoutGroundVeg - LWRabs_T.LWRabsGroundVeg;
LWREB_T.LWREBTree				=	LWRin_T.LWRinTree - LWRout_T.LWRoutTree - LWRabs_T.LWRabsTree;
LWREB_T.LWREBWallSun			=	LWRin_T.LWRinWallSun-LWRout_T.LWRoutWallSun - LWRabs_T.LWRabsWallSun;
LWREB_T.LWREBWallShade			=	LWRin_T.LWRinWallShade-LWRout_T.LWRoutWallShade - LWRabs_T.LWRabsWallShade;
LWREB_T.LWREBTotalGround		=	LWRin_T.LWRinTotalGround-LWRout_T.LWRoutTotalGround - LWRabs_T.LWRabsTotalGround;
LWREB_T.LWREBTotalCanyon		=	LWRin_T.LWRinTotalCanyon-LWRout_T.LWRoutTotalCanyon - LWRabs_T.LWRabsTotalCanyon;

if abs(LWREB_T.LWREBGroundImp)>=10^-6
	disp('LWREB_T.LWREBGroundImp is not 0. Please check LWRabsorbedWithTrees.m')
end
if abs(LWREB_T.LWREBGroundBare)>=10^-6
	disp('LWREB_T.LWREBGroundBare is not 0. Please check LWRabsorbedWithTrees.m')
end
if abs(LWREB_T.LWREBGroundVeg)>=10^-6
	disp('LWREB_T.LWREBGroundVeg	 is not 0. Please check LWRabsorbedWithTrees.m')
end
if abs(LWREB_T.LWREBWallSun)>=10^-6
	disp('LWREB_T.LWREBWallSun is not 0. Please check LWRabsorbedWithTrees.m')
end
if abs(LWREB_T.LWREBWallShade)>=10^-6
	disp('LWREB_T.LWREBWallShade is not 0. Please check LWRabsorbedWithTrees.m')
end
if abs(LWREB_T.LWREBTotalGround)>=10^-6
	disp('LWREB_T.LWREBTotalGround is not 0. Please check LWRabsorbedWithTrees.m')
end
if abs(LWREB_T.LWREBTotalCanyon)>=10^-6
	disp('LWREB_T.LWREBTotalCanyon is not 0. Please check LWRabsorbedWithTrees.m')
end
if abs(LWREB_T.LWREBTree)>=10^-6
	disp('LWREB_T.LWREBTree is not 0. Please check LWRabsorbedWithTrees.m')
end
