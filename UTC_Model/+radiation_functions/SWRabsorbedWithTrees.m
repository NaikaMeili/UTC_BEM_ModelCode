function[SWRin_T,SWRout_T,SWRabs_T,SWRabsDir_T,SWRabsDiff_T,SWREB_T,albedo_canyon]...
         =SWRabsorbedWithTrees(h_can,w_can,h_tree,r_tree,d_tree,...
         fgveg,fgbare,fgimp,awraw,agveg,agbare,agimp,at,...
         LAIt,SWR_dir,SWR_diff,theta_Z,theta_n,ViewFactor,ParVegTree,ParWindows,BEM_on)
	 
     
% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SWRin_T			=	total incoming shortwave radiation to surface i [W/m^2]
% SWRout_T			=	total outgoing shortwave radiation of surface i [W/m^2]
% SWRabs_T			=	total absorbed shortwave radiation of surface i [W/m^2]
% SWRabsDir_T		=	total absorbed direct shortwave radiation of surface i [W/m^2]
% SWRabsDiff_T		=	total absorbed diffuse shortwace radiation of surface i [W/m^2]
% SWREB_T			=	Surface energy balance of each surface i [W/m^2]

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h_can					=	building height [-]
% w_can					=	ground width [-]
% d_tree				=	location of trees in the canyon, tree-wall distance  [-]
% h_tree				=	height of trees, vertical level at the crown center  [-]
% r_tree				=	size of the tree crown, crown radius  [-]
% fgveg					=	Partitioning ground vegetation [-]
% fgbare				=	Partitioning ground bare [-]
% fgimp					=	Partitioning ground impervious [-]
% at					=	tree surface albedo  [-]
% aw					=	Wall surface albedo  [-]
% agveg					=	ground vegetation albedo  [-]
% agbare				=	ground bare albedo  [-]
% agimp					=	ground impervious albedo  [-]
% LAIt					=	leaf area index of the tree
% SWR_dir				=	direct shortwave radiation W/m^2 of horizontal surfaces [W/m^2]
% SWR_diff				=	diffuce shortwave radiation W per m^2 of horizontal surface [W/m^2]
% theta_Z				=	solar zenith angle [rad]
% theta_n				=	difference between solar azimuth angle and canyon orientation [rad]

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

% load shortwave radiation
[SWRdir_ground,SWRdir_wallsun,SWRdir_wallshade,SWRdir_tree]=...
    radiation_functions.DirectSWRSurfaces(h_can,w_can,d_tree,h_tree,r_tree,theta_Z,theta_n,SWR_dir,LAIt,1,ParVegTree);

% Balance direct shortwave radiation in
EB_SWRdir		=	SWR_dir-(SWRdir_ground*A_g/A_g+SWRdir_wallsun*A_w/A_g+SWRdir_wallshade*A_w/A_g+SWRdir_tree*A_t/A_g);
EB_SWRdiff		=	SWR_diff-(F_sg_T*SWR_diff+F_sw_T*SWR_diff+F_sw_T*SWR_diff+F_st_T*SWR_diff);

if abs(EB_SWRdir)>=10^-6
	disp('EB_SWRdir is not 0. Please check SWRabsorbedWithTrees.m')
end

if abs(EB_SWRdiff)>=10^-6
	disp('EB_SWRdiff is not 0. Please check SWRabsorbedWithTrees.m')
end

% Check if view factors add up to 1
SVF(1)			=	F_gs_T+F_gt_T+2*F_gw_T;
SVF(2)			=	F_ww_T+F_wt_T+F_wg_T+F_ws_T;
SVF(3)			=	F_sg_T+2*F_sw_T+F_st_T;
SVF(4)			=	F_ts_T+2*F_tw_T+F_tt_T+F_tg_T;

SVF2(1)			=	F_gs_T+2*F_ws_T*A_w+F_ts_T*A_t;
SVF2(2)			=	F_sg_T+2*F_wg_T*A_w+F_tg_T*A_t;
SVF2(3)			=	F_ww_T+F_sw_T*A_g/A_w+F_gw_T*A_g/A_w+F_tw_T*A_t/A_w;
SVF2(4)			=	F_gt_T*A_g/A_t+2*F_wt_T*A_w/A_t+F_tt_T;

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
% Trees
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cimp	=	fgimp>0;
Cbare	=	fgbare>0;
Cveg	=	fgveg>0;

if BEM_on==1
    aw = awraw.*(1-ParWindows.GlazingRatio) + ParWindows.SolarAlbedo.*ParWindows.GlazingRatio;
else
    aw = awraw;
end

% Albedos
ai=[agveg;agbare;agimp;aw;aw;at;0];

% View factor matrix to solve for infinite reflections equation
% Omega_i = Tij*Bj
Tij	=	[1,0,0, -agveg*F_gw_T*Cveg, -agveg*F_gw_T*Cveg, -agveg*F_gt_T*Cveg, -agveg*F_gs_T*Cveg;...
		0,1,0, -agbare*F_gw_T*Cbare, -agbare*F_gw_T*Cbare, -agbare*F_gt_T*Cbare, -agbare*F_gs_T*Cbare;...
		0,0,1, -agimp*F_gw_T*Cimp, -agimp*F_gw_T*Cimp, -agimp*F_gt_T*Cimp, -agimp*F_gs_T*Cimp;...
		-aw*F_wg_T*fgveg*Cveg,-aw*F_wg_T*fgbare*Cbare,-aw*F_wg_T*fgimp*Cimp, 1, -aw*F_ww_T, -aw*F_wt_T, -aw*F_ws_T;...
		-aw*F_wg_T*fgveg*Cveg,-aw*F_wg_T*fgbare*Cbare,-aw*F_wg_T*fgimp*Cimp, -aw*F_ww_T, 1, -aw*F_wt_T, -aw*F_ws_T;...
		-at*F_tg_T*fgveg*Cveg,-at*F_tg_T*fgbare*Cbare,-at*F_tg_T*fgimp*Cimp, -at*F_tw_T, -at*F_tw_T, 1-at*F_tt_T, -at*F_ts_T;...
		0, 0, 0, 0, 0, 0, 1];

% Incoming shortwave radiation from sky
Omega_i	=	[agveg*SWRdir_ground*Cveg;...
			agbare*SWRdir_ground*Cbare;...
			agimp*SWRdir_ground*Cimp;...
			aw*SWRdir_wallsun;...
			aw*0;...
			at*SWRdir_tree;...
			SWR_diff];
		
% How to solve the set of equations
% The outgoing and emitted radiation should be the same
% B_i				=	[Bveg; Bbare; Bimp; Bwall; Bwall; Bsky];
% Omega_i			=	Tij*B_i
% Tij^-1*Omega_i	=	Tij^-1*Tij*B_i
% Tij^-1*Omega_i	=	B_i

% Outgoing radiation per surface
B_i		=	Tij^-1*Omega_i;	% Outgoing radiation [W/m^2] per m^2 surface area

if B_i(7,1)~=SWR_diff
	disp('Incoming lonwave radiation and emitted longwave radiation from the sky after the matrix inversion are not equal')
end

% Incoming shortwave radiation at each surface A_i
Tij2=	[0, 0, 0, F_gw_T*Cveg, F_gw_T*Cveg, F_gt_T*Cveg, F_gs_T*Cveg;...
		0, 0, 0, F_gw_T*Cbare, F_gw_T*Cbare, F_gt_T*Cbare, F_gs_T*Cbare;...
		0, 0, 0, F_gw_T*Cimp, F_gw_T*Cimp, F_gt_T*Cimp, F_gs_T*Cimp;...
		F_wg_T*fgveg*Cveg, F_wg_T*fgbare*Cbare, F_wg_T*fgimp*Cimp, 0, F_ww_T, F_wt_T, F_ws_T;...
		F_wg_T*fgveg*Cveg, F_wg_T*fgbare*Cbare, F_wg_T*fgimp*Cimp, F_ww_T, 0, F_wt_T, F_ws_T;...
		F_tg_T*fgveg*Cveg, F_tg_T*fgbare*Cbare, F_tg_T*fgimp*Cimp, F_tw_T, F_tw_T, F_tt_T, F_ts_T;...
		0, 0, 0, 0, 0, 0, 0];
	
SWRdir_i=	[SWRdir_ground*Cveg;...
			SWRdir_ground*Cbare;...
			SWRdir_ground*Cimp;...
			SWRdir_wallsun;...
			0;...
			SWRdir_tree;...
			0];

A_i1	=	Tij2*B_i+SWRdir_i;	% Incoming radiation [W/m^2] per m^2 surface area

A_i			=	B_i./ai;		% Incoming radiation [W/m^2] per m^2 surface area
A_i(ai==0)	=	A_i1(ai==0);
A_i(7)		=	0;				% Assumption: The sky has a fixed emission of LWR. Hence, Qnet is 0.

% Absorbed shortwave radiation at ech surface Qnet_i
Qnet_i		=	A_i-B_i;

% Assignment
SWRout_i		=	B_i;	% Outgoing radiation [W/m^2] per m^2 surface area
SWRin_i			=	A_i;	% Incoming radiation [W/m^2] per m^2 surface area
SWRnet_i		=	Qnet_i;	% Net absorbed radiation [W/m^2] per m^2 surface area

% % Energy balance
SWRin_atm			=	SWR_dir+SWR_diff;

TotalSWRSurface_in	=	SWRin_i(1)*fgveg*A_g/A_g + SWRin_i(2)*fgbare*A_g/A_g + SWRin_i(3)*fgimp*A_g/A_g + ...
						SWRin_i(4)*A_w/A_g + SWRin_i(5)*A_w/A_g + SWRin_i(6)*A_t/A_g;

TotalSWRSurface_abs	=	SWRnet_i(1)*fgveg*A_g/A_g + SWRnet_i(2)*fgbare*A_g/A_g + SWRnet_i(3)*fgimp*A_g/A_g + ...
						SWRnet_i(4)*A_w/A_g + SWRnet_i(5)*A_w/A_g + SWRnet_i(6)*A_t/A_g;
      
TotalSWRSurface_out	=	SWRout_i(1)*fgveg*A_g/A_s+SWRout_i(2)*fgbare*A_g/A_s+SWRout_i(3)*fgimp*A_g/A_s+...
						SWRout_i(4)*A_w/A_s+SWRout_i(5)*A_w/A_s+SWRout_i(6)*A_t/A_s;
	  
TotalSWRref_to_atm	=	SWRout_i(1)*F_sg_T*fgveg + SWRout_i(2)*F_sg_T*fgbare + SWRout_i(3)*F_sg_T*fgimp + ...
						SWRout_i(4)*F_sw_T + SWRout_i(5)*F_sw_T + SWRout_i(6)*F_st_T;
					
TotalSWRref_to_atm2	=	SWRout_i(1)*F_gs_T*fgveg + SWRout_i(2)*F_gs_T*fgbare + SWRout_i(3)*F_gs_T*fgimp + ...
						SWRout_i(4)*F_ws_T + SWRout_i(5)*F_ws_T + SWRout_i(6)*F_ts_T;

EBSurface			=	TotalSWRSurface_in - TotalSWRSurface_abs - TotalSWRSurface_out;
EBCanyon			=	SWRin_atm - TotalSWRSurface_abs - TotalSWRref_to_atm;

albedo_canyon		=	TotalSWRref_to_atm/SWRin_atm;

% Energy balance
if abs(EBSurface)>=10^-6
	disp('EBSurface is not 0. Please check SWRabsorbedWithTrees.m')
end
if abs(EBCanyon)>=10^-6
	disp('EBCanyon is not 0. Please check SWRabsorbedWithTrees.m')
end

if abs(EBCanyon)>=10^-6
	msgbox('EBCanyon is not 0, calculate the appropriate view factors', 'Error','error');
	error('EBCanyon is not 0, calculate the appropriate view factors')	
end	

% Sequence in Vectors :
% Vegetated ground
% Bare ground
% Impervious ground
% Sunlit wall
% Shaded wall
% Trees
% Shortwave radiation by each surface per m^2 surface area
% Incoming shortwave radiation
SWRin_T							=	[];			
SWRin_T.SWRinGroundImp			=	SWRin_i(3)*Cimp;
SWRin_T.SWRinGroundBare			=	SWRin_i(2)*Cbare;
SWRin_T.SWRinGroundVeg			=	SWRin_i(1)*Cveg;
SWRin_T.SWRinTree				=	SWRin_i(6);
SWRin_T.SWRinWallSun			=	SWRin_i(4);
SWRin_T.SWRinWallShade			=	SWRin_i(5);
SWRin_T.SWRinTotalGround		=	fgveg*SWRin_i(1)+fgbare*SWRin_i(2)+fgimp*SWRin_i(3);
SWRin_T.SWRinTotalCanyon		=	SWRin_i(1)*fgveg*A_g/A_g + SWRin_i(2)*fgbare*A_g/A_g + SWRin_i(3)*fgimp*A_g/A_g + ...
									SWRin_i(4)*A_w/A_g + SWRin_i(5)*A_w/A_g + SWRin_i(6)*A_t/A_g;
								
% Outgoing shortwave radiation
SWRout_T						=	[];			
SWRout_T.SWRoutGroundImp		=	SWRout_i(3)*Cimp;
SWRout_T.SWRoutGroundBare		=	SWRout_i(2)*Cbare;
SWRout_T.SWRoutGroundVeg		=	SWRout_i(1)*Cveg;
SWRout_T.SWRoutTree				=	SWRout_i(6);
SWRout_T.SWRoutWallSun			=	SWRout_i(4);
SWRout_T.SWRoutWallShade		=	SWRout_i(5);
SWRout_T.SWRoutTotalGround		=	fgveg*SWRout_i(1)+fgbare*SWRout_i(2)+fgimp*SWRout_i(3);
SWRout_T.SWRoutTotalCanyon		=	SWRout_i(1)*fgveg*A_g/A_g + SWRout_i(2)*fgbare*A_g/A_g + SWRout_i(3)*fgimp*A_g/A_g + ...
									SWRout_i(4)*A_w/A_g + SWRout_i(5)*A_w/A_g + SWRout_i(6)*A_t/A_g;
								
% Absorbed shortwave radiation
SWRabs_T						=	[];			
SWRabs_T.SWRabsGroundImp		=	SWRnet_i(3)*Cimp;
SWRabs_T.SWRabsGroundBare		=	SWRnet_i(2)*Cbare;
SWRabs_T.SWRabsGroundVeg		=	SWRnet_i(1)*Cveg;
SWRabs_T.SWRabsTree				=	SWRnet_i(6);
SWRabs_T.SWRabsWallSun			=	SWRnet_i(4);
SWRabs_T.SWRabsWallShade		=	SWRnet_i(5);
SWRabs_T.SWRabsTotalGround		=	fgveg*SWRnet_i(1)+fgbare*SWRnet_i(2)+fgimp*SWRnet_i(3);
SWRabs_T.SWRabsTotalCanyon		=	SWRnet_i(1)*fgveg*A_g/A_g + SWRnet_i(2)*fgbare*A_g/A_g + SWRnet_i(3)*fgimp*A_g/A_g + ...
									SWRnet_i(4)*A_w/A_g + SWRnet_i(5)*A_w/A_g + SWRnet_i(6)*A_t/A_g;
								
% Direct absorbed shortwave radiation
SWRabsDir_T						=	[];
SWRabsDir_T.SWRabsGroundImp		=	(1-agimp)*SWRdir_ground*Cimp;
SWRabsDir_T.SWRabsGroundBare	=	(1-agbare)*SWRdir_ground*Cbare;
SWRabsDir_T.SWRabsGroundVeg		=	(1-agveg)*SWRdir_ground*Cveg;
SWRabsDir_T.SWRabsTree			=	(1-at)*SWRdir_tree;
SWRabsDir_T.SWRabsWallSun		=	(1-aw)*SWRdir_wallsun;
SWRabsDir_T.SWRabsWallShade		=	(1-aw)*SWRdir_wallshade;
SWRabsDir_T.SWRabsTotalGround	=	fgveg*(1-agveg)*SWRdir_ground+fgbare*(1-agbare)*SWRdir_ground+fgimp*(1-agimp)*SWRdir_ground;
SWRabsDir_T.SWRabsTotalCanyon	=	fgveg*(1-agveg)*SWRdir_ground*A_g/A_g+fgbare*(1-agbare)*SWRdir_ground*A_g/A_g+fgimp*(1-agimp)*SWRdir_ground*A_g/A_g + ...
									(1-aw)*SWRdir_wallsun*A_w/A_g + (1-aw)*SWRdir_wallshade*A_w/A_g + (1-at)*SWRdir_tree*A_t/A_g;

% Diffuse absorbed shortwave radiation
SWRabsDiff_T					=	[];			% Absorbed shortwave radiation
SWRabsDiff_T.SWRabsGroundImp	=	(SWRabs_T.SWRabsGroundImp-SWRabsDir_T.SWRabsGroundImp)*Cimp;
SWRabsDiff_T.SWRabsGroundBare	=	(SWRabs_T.SWRabsGroundBare-SWRabsDir_T.SWRabsGroundBare)*Cbare;
SWRabsDiff_T.SWRabsGroundVeg	=	(SWRabs_T.SWRabsGroundVeg-SWRabsDir_T.SWRabsGroundVeg)*Cveg;
SWRabsDiff_T.SWRabsTree			=	SWRabs_T.SWRabsTree-SWRabsDir_T.SWRabsTree;
SWRabsDiff_T.SWRabsWallSun		=	SWRabs_T.SWRabsWallSun-SWRabsDir_T.SWRabsWallSun;
SWRabsDiff_T.SWRabsWallShade	=	SWRabs_T.SWRabsWallShade-SWRabsDir_T.SWRabsWallShade;
SWRabsDiff_T.SWRabsTotalGround	=	SWRabs_T.SWRabsTotalGround-SWRabsDir_T.SWRabsTotalGround;
SWRabsDiff_T.SWRabsTotalCanyon	=	SWRabs_T.SWRabsTotalCanyon-SWRabsDir_T.SWRabsTotalCanyon;

% Energy Balance of shortwave radiation								
SWREB_T							=	[];			
SWREB_T.SWREBGroundImp			=	SWRin_T.SWRinGroundImp - SWRout_T.SWRoutGroundImp - SWRabs_T.SWRabsGroundImp;
SWREB_T.SWREBGroundBare			=	SWRin_T.SWRinGroundBare - SWRout_T.SWRoutGroundBare - SWRabs_T.SWRabsGroundBare;
SWREB_T.SWREBGroundVeg			=	SWRin_T.SWRinGroundVeg - SWRout_T.SWRoutGroundVeg - SWRabs_T.SWRabsGroundVeg;
SWREB_T.SWREBTree				=	SWRin_T.SWRinTree - SWRout_T.SWRoutTree - SWRabs_T.SWRabsTree;
SWREB_T.SWREBWallSun			=	SWRin_T.SWRinWallSun-SWRout_T.SWRoutWallSun - SWRabs_T.SWRabsWallSun;
SWREB_T.SWREBWallShade			=	SWRin_T.SWRinWallShade-SWRout_T.SWRoutWallShade - SWRabs_T.SWRabsWallShade;
SWREB_T.SWREBTotalGround		=	SWRin_T.SWRinTotalGround-SWRout_T.SWRoutTotalGround - SWRabs_T.SWRabsTotalGround;
SWREB_T.SWREBTotalCanyon		=	SWRin_T.SWRinTotalCanyon-SWRout_T.SWRoutTotalCanyon - SWRabs_T.SWRabsTotalCanyon;


if abs(SWREB_T.SWREBGroundImp)>=10^-6
	disp('SWREB_T.SWREBGroundImp is not 0. Please check SWRabsorbedNoTrees.m')
elseif abs(SWREB_T.SWREBGroundBare)>=10^-6
	disp('SWREB_T.SWREBGroundBare is not 0. Please check SWRabsorbedNoTrees.m')
elseif abs(SWREB_T.SWREBGroundVeg)>=10^-6
	disp('SWREB_T.SWREBGroundVeg	 is not 0. Please check SWRabsorbedNoTrees.m')
elseif abs(SWREB_T.SWREBWallSun)>=10^-6
	disp('SWREB_T.SWREBWallSun is not 0. Please check SWRabsorbedNoTrees.m')
elseif abs(SWREB_T.SWREBWallShade)>=10^-6
	disp('SWREB_T.SWREBWallShade is not 0. Please check SWRabsorbedNoTrees.m')
elseif abs(SWREB_T.SWREBTotalGround)>=10^-6
	disp('SWREB_T.SWREBTotalGround is not 0. Please check SWRabsorbedNoTrees.m')
elseif abs(SWREB_T.SWREBTotalCanyon)>=10^-6
	disp('SWREB_T.SWREBTotalCanyon is not 0. Please check SWRabsorbedNoTrees.m')
elseif abs(SWREB_T.SWREBTree)>=10^-6
	disp('SWREB_T.SWREBTree is not 0. Please check SWRabsorbedNoTrees.m')
end



