function[LWRinB,LWRoutB,LWRabsB]=LWRabsBuildingHalf(Tceiling,Tinwall,Tintmass,Tground,...
    A_c,A_g,A_h,F_gc,F_gw,F_ww,F_wg,F_wc,F_cg,F_cw,ec,eg,ew,em)


%--------------------------------------------
bolzm	=	5.67*10^(-8);	% Stefan-Boltzmann constant [W*m^-2*K-4]

% Calculations for infinite reflections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sequence in Vectors :
% Ceiling
% Wall
% Internal mass (internal wall)
% Ground
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% View factor matrix to solve for infinite reflections equation
% Omega_i = Tij*Bj
% B_i		=	[Bceiling; Bwallsun; Bwallshade; Bground];
Tij	=	[1, -(1-ec)*F_cw, -(1-ec)*F_cw, -(1-ec)*F_cg;...
         -(1-ew)*F_wc, 1, -(1-ew)*F_ww, -(1-ew)*F_wg;...
         -(1-em)*F_wc, -(1-em)*F_ww, 1, -(1-em)*F_wg;...
         -(1-eg)*F_gc, -(1-eg)*F_gw, -(1-eg)*F_gw, 1];
	
% Emitted radiation per surface
Omega_i	=	[(ec*bolzm*(Tceiling)^4);...
			(ew*bolzm*(Tinwall)^4);...
			(em*bolzm*(Tintmass)^4);...
			(eg*bolzm*(Tground)^4)];
		
% How to solve the set of equations
% The outgoing and emitted radiation should be the same
% B_i				=	[Bceiling; Bwallsun; Bwallshade; Bground];
% Omega_i			=	Tij*B_i
% Tij^-1*Omega_i	=	Tij^-1*Tij*B_i
% Tij^-1*Omega_i	=	B_i

% Outgoing radiation per surface
B_i		=	Tij^-1*Omega_i;	% Outgoing radiation [W/m^2] per m^2 surface area

% Incoming longwave radiation at each surface A_i
Tij2	=	[0, F_cw, F_cw, F_cg;...
         F_wc, 0, F_ww, F_wg;...
         F_wc, F_ww, 0, F_wg;...
         F_gc, F_gw, F_gw, 0];

A_i		=	Tij2*B_i;
e_i		=	[ec; ew; em; eg];
A_i2	=	(B_i-Omega_i)./(1-e_i);
Qnet_i2	=	A_i-B_i;


% Absorbed longwave radiation (Harman et al 2004)
Qnet_i			=	(e_i.*B_i - Omega_i)./(1-e_i);
Qnet_i(e_i==1)	=	A_i(e_i==1) - Omega_i(e_i==1);


% Difference between emitted and outgoing
% Delta_out		=	B_i-Omega_i;

% Assignment
LWRout_i		=	B_i;	% Outgoing radiation [W/m^2] per m^2 surface area
LWRemit_i		=	Omega_i;% Emitted radiation [W/m^2] per m^2 surface area
LWRin_i			=	A_i;	% Incoming radiation [W/m^2] per m^2 surface area
LWRnet_i		=	Qnet_i;	% Net absorbed radiation [W/m^2] per m^2 surface area

% Energy balance
TotalLWRSurface_in	=	LWRin_i(1)*A_c/A_g + LWRin_i(2)*A_h/A_g + LWRin_i(3)*A_h/A_g + LWRin_i(4)*A_g/A_g;

TotalLWRSurface_abs	=	LWRnet_i(1)*A_c/A_g + LWRnet_i(2)*A_h/A_g + LWRnet_i(3)*A_h/A_g + LWRnet_i(4)*A_g/A_g;
      
TotalLWRSurface_out	=	LWRout_i(1)*A_c/A_g + LWRout_i(2)*A_h/A_g + LWRout_i(3)*A_h/A_g + LWRout_i(4)*A_g/A_g;

EBSurface			=	TotalLWRSurface_in - TotalLWRSurface_abs - TotalLWRSurface_out;

% Energy balance
if abs(EBSurface)>=10^-6
	disp('EBSurface is not 0. Please check LWRabsIndoors.m')
end


% Longwave radiation by each surface per m^2 surface area
% Incoming longwave radiation
LWRinB                  =	[];			
LWRinB.LWRinCeiling     =	LWRin_i(1);
LWRinB.LWRinWall        =	LWRin_i(2);
LWRinB.LWRinInternalMass=	LWRin_i(3);
LWRinB.LWRinGround      =	LWRin_i(4);

								
% Outgoing longwave radiation
LWRoutB						=	[];			
LWRoutB.LWRoutCeiling		=	LWRout_i(1);
LWRoutB.LWRoutWall		    =	LWRout_i(2);
LWRoutB.LWRoutInternalMass  =	LWRout_i(3);
LWRoutB.LWRoutGround		=	LWRout_i(4);


% Absorbed longwave radiation
LWRabsB						=	[];			
LWRabsB.LWRabsCeiling		=	LWRnet_i(1);
LWRabsB.LWRabsWall		    =	LWRnet_i(2);
LWRabsB.LWRabsInternalMass  =	LWRnet_i(3);
LWRabsB.LWRabsGround        =	LWRnet_i(4);


% Energy Balance of longwave radiation								
LWREBB				        =	[];			
LWREBB.LWREBCeiling         =	LWRinB.LWRinCeiling - LWRoutB.LWRoutCeiling - LWRabsB.LWRabsCeiling;
LWREBB.LWREBWall		    =	LWRinB.LWRinWall - LWRoutB.LWRoutWall - LWRabsB.LWRabsWall;
LWREBB.LWREBInternalMass    =	LWRinB.LWRinInternalMass - LWRoutB.LWRoutInternalMass - LWRabsB.LWRabsInternalMass;
LWREBB.LWREBGround          =	LWRinB.LWRinGround - LWRoutB.LWRoutGround - LWRabsB.LWRabsGround;



if abs(LWREBB.LWREBCeiling)>=10^-6 
	disp('LWREBB.LWREBCeiling is not 0. Please check LWRabsIndoors.m')
end
if abs(LWREBB.LWREBWall)>=10^-6
	disp('LWREBB.LWREBWall is not 0. Please check LWRabsIndoors.m')
end
if abs(LWREBB.LWREBInternalMass )>=10^-6
	disp('LWREBB.LWREBInternalMass	 is not 0. Please check LWRabsIndoors.m')
end
if abs(LWREBB.LWREBGround)>=10^-6
	disp('LWREBB.LWREBGround is not 0. Please check LWRabsIndoors.m')
end

