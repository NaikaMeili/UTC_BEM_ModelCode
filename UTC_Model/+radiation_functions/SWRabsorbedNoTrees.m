function[SWRin_nT,SWRout_nT,SWRabs_nT,SWRabsDir_nT,SWRabsDiff_nT,SWREB_nT,albedo_canyon]...
         =SWRabsorbedNoTrees(h_can,w_can,fgveg,fgbare,fgimp,awraw,agveg,agbare,agimp,...
         SWR_dir,SWR_diff,theta_Z,theta_n,ViewFactor,ParVegTree,ParWindows,BEM_on)
 
% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SWRin_nT			=	total incoming shortwave radiation to surface i [W/m^2]
% SWRout_nT			=	total outgoing shortwave radiation of surface i [W/m^2]
% SWRabs_nT			=	total absorbed shortwave radiation of surface i [W/m^2]
% SWRabsDir_nT		=	total absorbed direct shortwave radiation of surface i [W/m^2]
% SWRabsDiff_nT		=	total absorbed diffuse shortwace radiation of surface i [W/m^2]
% SWREB_nT			=	Surface energy balance of each surface i [W/m^2]

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h_can						=	building height [-]
% w_can						=	ground width [-]
% fgveg						=	Partitioning ground vegetation [-]
% fgbare					=	Partitioning ground bare [-]
% fgimp						=	Partitioning ground impervious [-]
% aw						=	Wall surface albedo [-]
% agveg						=	ground vegetation albedo [-]
% agbare					=	ground bare albedo [-]
% agimp						=	ground impervious albedo [-]
% SWR_dir					=	direct shortwave radiation W/m^2 of horizontal surfaces  [W/m^2]
% SWR_diff					=	diffuce shortwave radiation W per m^2 of horizontal surface  [W/m^2]
% theta_Z					=	solar zenith angle [rad]
% theta_n					=	difference between solar azimuth angle and canyon orientation [rad]

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

% Check if view factors add up to 1
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

% load shortwave radiation
[SWRdir_ground,SWRdir_wallsun,SWRdir_wallshade,~]=...
    radiation_functions.DirectSWRSurfaces(h_can,w_can,NaN,NaN,NaN,theta_Z,theta_n,SWR_dir,NaN,0,ParVegTree);

% Balance direct shortwave radiation in
EB_SWRdir		=	SWR_dir-(SWRdir_ground*A_g/A_g+SWRdir_wallsun*A_w/A_g+SWRdir_wallshade*A_w/A_g);
EB_SWRdiff		=	SWR_diff-(F_sg_nT*SWR_diff+F_sw_nT*SWR_diff+F_sw_nT*SWR_diff);

if abs(EB_SWRdir)>=10^-6
	disp('EB_SWRdir is not 0. Please check SWRabsorbedNoTrees.m')
end
if abs(EB_SWRdiff)>=10^-6
	disp('EB_SWRdiff is not 0. Please check SWRabsorbedNoTrees.m')
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
Cimp	=	fgimp>0;
Cbare	=	fgbare>0;
Cveg	=	fgveg>0;

if BEM_on==1
    aw = awraw.*(1-ParWindows.GlazingRatio) + ParWindows.SolarAlbedo.*ParWindows.GlazingRatio;
else
    aw = awraw;
end

% Albedos
ai	=	[agveg;agbare;agimp;aw;aw;0];

% View factor matrix to solve for infinite reflections equation
% Omega_i = Tij*Bj
% B_i		=	[Bveg; Bbare; Bimp; Bwall; Bwall; Bsky];
Tij	=	[1,0,0, -agveg*F_gw_nT*Cveg, -agveg*F_gw_nT*Cveg, -agveg*F_gs_nT*Cveg;...
		0,1,0, -agbare*F_gw_nT*Cbare, -agbare*F_gw_nT*Cbare, -agbare*F_gs_nT*Cbare;...
		0,0,1, -agimp*F_gw_nT*Cimp, -agimp*F_gw_nT*Cimp, -agimp*F_gs_nT*Cimp;...
		-aw*F_wg_nT*fgveg*Cveg,-aw*F_wg_nT*fgbare*Cbare,-aw*F_wg_nT*fgimp*Cimp, 1, -aw*F_ww_nT, -aw*F_ws_nT;...
		-aw*F_wg_nT*fgveg*Cveg,-aw*F_wg_nT*fgbare*Cbare,-aw*F_wg_nT*fgimp*Cimp, -aw*F_ww_nT, 1, -aw*F_ws_nT;
		0, 0, 0, 0, 0, 1];

% B_i		=	[Bveg; Bbare; Bimp; Bwall; Bwall; Bsky];
Tij	=	[1,0,0, -agveg*F_gw_nT*Cveg, -agveg*F_gw_nT*Cveg, -agveg*F_gs_nT*Cveg;...
		0,1,0, -agbare*F_gw_nT*Cbare, -agbare*F_gw_nT*Cbare, -agbare*F_gs_nT*Cbare;...
		0,0,1, -agimp*F_gw_nT*Cimp, -agimp*F_gw_nT*Cimp, -agimp*F_gs_nT*Cimp;...
		-aw*F_wg_nT*fgveg*Cveg,-aw*F_wg_nT*fgbare*Cbare,-aw*F_wg_nT*fgimp*Cimp, 1, -aw*F_ww_nT, -aw*F_ws_nT;...
		-aw*F_wg_nT*fgveg*Cveg,-aw*F_wg_nT*fgbare*Cbare,-aw*F_wg_nT*fgimp*Cimp, -aw*F_ww_nT, 1, -aw*F_ws_nT;
		0, 0, 0, 0, 0, 1];
	
% Incoming shortwave radiation from sky
Omega_i	=	[agveg*SWRdir_ground*Cveg;...
			agbare*SWRdir_ground*Cbare;...
			agimp*SWRdir_ground*Cimp;...
			aw*SWRdir_wallsun;...
			aw*0;...
			SWR_diff];
	
% How to solve the set of equations
% The outgoing and emitted radiation should be the same
% B_i				=	[Bveg; Bbare; Bimp; Bwall; Bwall; Bsky];
% Omega_i			=	Tij*B_i
% Tij^-1*Omega_i	=	Tij^-1*Tij*B_i
% Tij^-1*Omega_i	=	B_i

% Outgoing radiation per surface
B_i		=	Tij^-1*Omega_i;	% Outgoing radiation [W/m^2] per m^2 surface area

if B_i(6,1)~=SWR_diff
	disp('Incoming lonwave radiation and emitted longwave radiation from the sky after the matrix inversion are not equal')
end

% Incoming shortwave radiation at each surface A_i
Tij2	=	[0, 0, 0, F_gw_nT*Cveg, F_gw_nT*Cveg, F_gs_nT*Cveg;...
			0, 0, 0, F_gw_nT*Cbare, F_gw_nT*Cbare, F_gs_nT*Cbare;...
			0, 0, 0, F_gw_nT*Cimp, F_gw_nT*Cimp, F_gs_nT*Cimp;...
			F_wg_nT*fgveg*Cveg, F_wg_nT*fgbare*Cbare, F_wg_nT*fgimp*Cimp, 0, F_ww_nT, F_ws_nT;...
			F_wg_nT*fgveg*Cveg, F_wg_nT*fgbare*Cbare, F_wg_nT*fgimp*Cimp, F_ww_nT, 0, F_ws_nT;...
			0, 0, 0, 0, 0, 0];
		
SWRdir_i=	[SWRdir_ground*Cveg;...
			SWRdir_ground*Cbare;...
			SWRdir_ground*Cimp;...
			SWRdir_wallsun;...
			0;...
			0];

A_i1	=	Tij2*B_i+SWRdir_i;	% Incoming radiation [W/m^2] per m^2 surface area

A_i			=	B_i./ai;		% Incoming radiation [W/m^2] per m^2 surface area
A_i(ai==0)	=	A_i1(ai==0);
A_i(6)		=	0;				% Assumption: The sky has a fixed emission of LWR. Hence, Qnet is 0.

% Absorbed shortwave radiation at ech surface Qnet_i
Qnet_i		=	A_i-B_i;

% Assignment
SWRout_i		=	B_i;	% Outgoing radiation [W/m^2] per m^2 surface area
SWRin_i			=	A_i;	% Incoming radiation [W/m^2] per m^2 surface area
SWRnet_i		=	Qnet_i;	% Net absorbed radiation [W/m^2] per m^2 surface area

% % Energy balance
SWRin_atm			=	SWR_dir+SWR_diff;

TotalSWRSurface_in	=	SWRin_i(1)*fgveg*A_g/A_g + SWRin_i(2)*fgbare*A_g/A_g + SWRin_i(3)*fgimp*A_g/A_g + ...
						SWRin_i(4)*A_w/A_g + SWRin_i(5)*A_w/A_g;

TotalSWRSurface_abs	=	SWRnet_i(1)*fgveg*A_g/A_g + SWRnet_i(2)*fgbare*A_g/A_g + SWRnet_i(3)*fgimp*A_g/A_g + ...
						SWRnet_i(4)*A_w/A_g + SWRnet_i(5)*A_w/A_g;
      
TotalSWRSurface_out	=	SWRout_i(1)*fgveg*A_g/A_s+SWRout_i(2)*fgbare*A_g/A_s+SWRout_i(3)*fgimp*A_g/A_s+...
						SWRout_i(4)*A_w/A_s+SWRout_i(5)*A_w/A_s;
	  
TotalSWRref_to_atm	=	SWRout_i(1)*F_sg_nT*fgveg + SWRout_i(2)*F_sg_nT*fgbare + SWRout_i(3)*F_sg_nT*fgimp + ...
						SWRout_i(4)*F_sw_nT + SWRout_i(5)*F_sw_nT;
					
EBSurface			=	TotalSWRSurface_in - TotalSWRSurface_abs - TotalSWRSurface_out;
EBCanyon			=	SWRin_atm - TotalSWRSurface_abs - TotalSWRref_to_atm;

albedo_canyon		=	TotalSWRref_to_atm/SWRin_atm;

% Energy balance
if abs(EBSurface)>=10^-6
	disp('EBSurface is not 0. Please check SWRabsorbedNoTrees.m')
end
if abs(EBCanyon)>=10^-6
	disp('EBCanyon is not 0. Please check SWRabsorbedNoTrees.m')
end

% Sequence in Vectors :
% Vegetated ground
% Bare ground
% Impervious ground
% Sunlit wall
% Shaded wall

% Shortwave radiation by each surface per m^2 surface area
% Incoming shortwave radiation
SWRin_nT						=	[];			
SWRin_nT.SWRinGroundImp			=	SWRin_i(3)*Cimp;
SWRin_nT.SWRinGroundBare		=	SWRin_i(2)*Cbare;
SWRin_nT.SWRinGroundVeg			=	SWRin_i(1)*Cveg;
SWRin_nT.SWRinTree				=	0;
SWRin_nT.SWRinWallSun			=	SWRin_i(4);
SWRin_nT.SWRinWallShade			=	SWRin_i(5);
SWRin_nT.SWRinTotalGround		=	fgveg*SWRin_i(1)+fgbare*SWRin_i(2)+fgimp*SWRin_i(3);
SWRin_nT.SWRinTotalCanyon		=	SWRin_i(1)*fgveg*A_g/A_g + SWRin_i(2)*fgbare*A_g/A_g + SWRin_i(3)*fgimp*A_g/A_g + ...
									SWRin_i(4)*A_w/A_g + SWRin_i(5)*A_w/A_g;
								
% Outgoing shortwave radiation
SWRout_nT						=	[];			
SWRout_nT.SWRoutGroundImp		=	SWRout_i(3)*Cimp;
SWRout_nT.SWRoutGroundBare		=	SWRout_i(2)*Cbare;
SWRout_nT.SWRoutGroundVeg		=	SWRout_i(1)*Cveg;
SWRout_nT.SWRoutTree			=	0;
SWRout_nT.SWRoutWallSun			=	SWRout_i(4);
SWRout_nT.SWRoutWallShade		=	SWRout_i(5);
SWRout_nT.SWRoutTotalGround		=	fgveg*SWRout_i(1)+fgbare*SWRout_i(2)+fgimp*SWRout_i(3);
SWRout_nT.SWRoutTotalCanyon		=	SWRout_i(1)*fgveg*A_g/A_g + SWRout_i(2)*fgbare*A_g/A_g + SWRout_i(3)*fgimp*A_g/A_g + ...
									SWRout_i(4)*A_w/A_g + SWRout_i(5)*A_w/A_g;
								
% Absorbed shortwave radiation
SWRabs_nT						=	[];			
SWRabs_nT.SWRabsGroundImp		=	SWRnet_i(3)*Cimp;
SWRabs_nT.SWRabsGroundBare		=	SWRnet_i(2)*Cbare;
SWRabs_nT.SWRabsGroundVeg		=	SWRnet_i(1)*Cveg;
SWRabs_nT.SWRabsTree			=	0;
SWRabs_nT.SWRabsWallSun			=	SWRnet_i(4);
SWRabs_nT.SWRabsWallShade		=	SWRnet_i(5);
SWRabs_nT.SWRabsTotalGround		=	fgveg*SWRnet_i(1)+fgbare*SWRnet_i(2)+fgimp*SWRnet_i(3);
SWRabs_nT.SWRabsTotalCanyon		=	SWRnet_i(1)*fgveg*A_g/A_g + SWRnet_i(2)*fgbare*A_g/A_g + SWRnet_i(3)*fgimp*A_g/A_g + ...
									SWRnet_i(4)*A_w/A_g + SWRnet_i(5)*A_w/A_g;
								
% Direct absorbed shortwave radiation
SWRabsDir_nT					=	[];
SWRabsDir_nT.SWRabsGroundImp	=	(1-agimp)*SWRdir_ground*Cimp;
SWRabsDir_nT.SWRabsGroundBare	=	(1-agbare)*SWRdir_ground*Cbare;
SWRabsDir_nT.SWRabsGroundVeg	=	(1-agveg)*SWRdir_ground*Cveg;
SWRabsDir_nT.SWRabsTree			=	0;
SWRabsDir_nT.SWRabsWallSun		=	(1-aw)*SWRdir_wallsun;
SWRabsDir_nT.SWRabsWallShade	=	(1-aw)*SWRdir_wallshade;
SWRabsDir_nT.SWRabsTotalGround	=	fgveg*(1-agveg)*SWRdir_ground+fgbare*(1-agbare)*SWRdir_ground+fgimp*(1-agimp)*SWRdir_ground;
SWRabsDir_nT.SWRabsTotalCanyon	=	fgveg*(1-agveg)*SWRdir_ground*A_g/A_g+fgbare*(1-agbare)*SWRdir_ground*A_g/A_g+fgimp*(1-agimp)*SWRdir_ground*A_g/A_g + ...
									(1-aw)*SWRdir_wallsun*A_w/A_g + (1-aw)*SWRdir_wallshade*A_w/A_g;

% Diffuse absorbed shortwave radiation
SWRabsDiff_nT					=	[];			% Absorbed shortwave radiation
SWRabsDiff_nT.SWRabsGroundImp	=	(SWRabs_nT.SWRabsGroundImp-SWRabsDir_nT.SWRabsGroundImp)*Cimp;
SWRabsDiff_nT.SWRabsGroundBare	=	(SWRabs_nT.SWRabsGroundBare-SWRabsDir_nT.SWRabsGroundBare)*Cbare;
SWRabsDiff_nT.SWRabsGroundVeg	=	(SWRabs_nT.SWRabsGroundVeg-SWRabsDir_nT.SWRabsGroundVeg)*Cveg;
SWRabsDiff_nT.SWRabsTree		=	0;
SWRabsDiff_nT.SWRabsWallSun		=	SWRabs_nT.SWRabsWallSun-SWRabsDir_nT.SWRabsWallSun;
SWRabsDiff_nT.SWRabsWallShade	=	SWRabs_nT.SWRabsWallShade-SWRabsDir_nT.SWRabsWallShade;
SWRabsDiff_nT.SWRabsTotalGround	=	SWRabs_nT.SWRabsTotalGround-SWRabsDir_nT.SWRabsTotalGround;
SWRabsDiff_nT.SWRabsTotalCanyon	=	SWRabs_nT.SWRabsTotalCanyon-SWRabsDir_nT.SWRabsTotalCanyon;

% Energy Balance of shortwave radiation								
SWREB_nT						=	[];			
SWREB_nT.SWREBGroundImp			=	SWRin_nT.SWRinGroundImp - SWRout_nT.SWRoutGroundImp - SWRabs_nT.SWRabsGroundImp;
SWREB_nT.SWREBGroundBare		=	SWRin_nT.SWRinGroundBare - SWRout_nT.SWRoutGroundBare - SWRabs_nT.SWRabsGroundBare;
SWREB_nT.SWREBGroundVeg			=	SWRin_nT.SWRinGroundVeg - SWRout_nT.SWRoutGroundVeg - SWRabs_nT.SWRabsGroundVeg;
SWREB_nT.SWREBTree				=	0;
SWREB_nT.SWREBWallSun			=	SWRin_nT.SWRinWallSun-SWRout_nT.SWRoutWallSun - SWRabs_nT.SWRabsWallSun;
SWREB_nT.SWREBWallShade			=	SWRin_nT.SWRinWallShade-SWRout_nT.SWRoutWallShade - SWRabs_nT.SWRabsWallShade;
SWREB_nT.SWREBTotalGround		=	SWRin_nT.SWRinTotalGround-SWRout_nT.SWRoutTotalGround - SWRabs_nT.SWRabsTotalGround;
SWREB_nT.SWREBTotalCanyon		=	SWRin_nT.SWRinTotalCanyon-SWRout_nT.SWRoutTotalCanyon - SWRabs_nT.SWRabsTotalCanyon;

if abs(SWREB_nT.SWREBGroundImp)>=10^-6
	disp('SWREB_nT.SWREBGroundImp is not 0. Please check SWRabsorbedNoTrees.m')
elseif abs(SWREB_nT.SWREBGroundBare)>=10^-6
	disp('SWREB_nT.SWREBGroundBare is not 0. Please check SWRabsorbedNoTrees.m')
elseif abs(SWREB_nT.SWREBGroundVeg	)>=10^-6
	disp('SWREB_nT.SWREBGroundVeg	 is not 0. Please check SWRabsorbedNoTrees.m')
elseif abs(SWREB_nT.SWREBWallSun)>=10^-6
	disp('SWREB_nT.SWREBWallSun is not 0. Please check SWRabsorbedNoTrees.m')
elseif abs(SWREB_nT.SWREBWallShade)>=10^-6
	disp('SWREB_nT.SWREBWallShade is not 0. Please check SWRabsorbedNoTrees.m')
elseif abs(SWREB_nT.SWREBTotalGround)>=10^-6
	disp('SWREB_nT.SWREBTotalGround is not 0. Please check SWRabsorbedNoTrees.m')
elseif abs(SWREB_nT.SWREBTotalCanyon)>=10^-6
	disp('SWREB_nT.SWREBTotalCanyon is not 0. Please check SWRabsorbedNoTrees.m')
end


