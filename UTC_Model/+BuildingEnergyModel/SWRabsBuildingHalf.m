function[SWRinB,SWRoutB,SWRabsB]=SWRabsBuildingHalf(A_c,A_g,A_h,SWRinW,...
    F_gc,F_gw,F_ww,F_wg,F_wc,F_cg,F_cw,abc,abw,abm,abg)



% Calculations for infinite reflections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sequence in Vectors :
% Ceiling
% Wall
% Internal mass (internal wall)
% Ground
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Albedos
ai	=	[abc;abw;abm;abg];

% View factor matrix to solve for infinite reflections equation
% Omega_i = Tij*Bj
% B_i		=	[Bceiling; Bwall; Binternalmass; Bground];

Tij	=	[1, -abc*F_cw, -abc*F_cw, -abc*F_cg;...
		 -abw*F_wc, 1, -abw*F_ww, -abw*F_wg;...
         -abm*F_wc, -abm*F_ww, 1, -abm*F_wc;...
         -abg*F_gc, -abg*F_gw, -abg*F_gw, 1];
	
% Incoming shortwave radiation from sky
Omega_i	=	[abc*F_cw*SWRinW;...
			0;...
			abm*F_ww*SWRinW;...
			abg*F_gw*SWRinW];
	
% How to solve the set of equations
% The outgoing and emitted radiation should be the same
% B_i				=	[Bceiling; Bwallsun; Bwallshade; Bground];
% Omega_i			=	Tij*B_i
% Tij^-1*Omega_i	=	Tij^-1*Tij*B_i
% Tij^-1*Omega_i	=	B_i

% Outgoing radiation per surface
B_i		=	Tij^-1*Omega_i;	% Outgoing radiation [W/m^2] per m^2 surface area


% Incoming shortwave radiation at each surface A_i
Tij2	=	[0, F_cw, F_cw, F_cg;...
		    F_wc, 0, F_ww, F_wg;...
            F_wc, F_ww, 0, F_wc;...
            F_gc, F_gw, F_gw, 0];
		
SWRdir_i=	[F_cw*SWRinW;...
			0;...
			F_ww*SWRinW;...
			F_gw*SWRinW];

A_i1	=	Tij2*B_i+SWRdir_i;	% Incoming radiation [W/m^2] per m^2 surface area

A_i			=	B_i./ai;		% Incoming radiation [W/m^2] per m^2 surface area
A_i(ai==0)	=	A_i1(ai==0);

% Absorbed shortwave radiation at ech surface Qnet_i
Qnet_i		=	A_i-B_i;

% Assignment
SWRout_i		=	B_i;	% Outgoing radiation [W/m^2] per m^2 surface area
SWRin_i			=	A_i;	% Incoming radiation [W/m^2] per m^2 surface area
SWRnet_i		=	Qnet_i;	% Net absorbed radiation [W/m^2] per m^2 surface area

% % Energy balance
TotalSWRSurface_in	=	SWRin_i(1)*A_c/A_g + SWRin_i(2)*A_h/A_g + SWRin_i(3)*A_h/A_g + SWRin_i(4)*A_g/A_g;

TotalSWRSurface_abs	=	SWRnet_i(1)*A_c/A_g + SWRnet_i(2)*A_h/A_g + SWRnet_i(3)*A_h/A_g + SWRnet_i(4)*A_g/A_g;
      
TotalSWRSurface_out	=	SWRout_i(1)*A_c/A_g+SWRout_i(2)*A_h/A_g+SWRout_i(3)*A_h/A_g + SWRout_i(4)*A_g/A_g;
					
EBSurface			=	TotalSWRSurface_in - TotalSWRSurface_abs - TotalSWRSurface_out;

EBindoors = A_h/A_g.*(SWRinW) - TotalSWRSurface_abs;

% Energy balance
if abs(EBindoors)>=10^-6
	disp('EBindoors is not 0. Please check SWRabsIndoors.m')
end


% Energy balance
if abs(EBSurface)>=10^-6
	disp('EBSurface is not 0. Please check SWRabsIndoors.m')
end



% Shortwave radiation by each surface per m^2 surface area
% Incoming shortwave radiation
SWRinB           =	[];			
SWRinB.SWRinCeiling     =	SWRin_i(1);
SWRinB.SWRinWall        =	SWRin_i(2);
SWRinB.SWRinInternalMass=	SWRin_i(3);
SWRinB.SWRinGround      =	SWRin_i(4);

								
% Outgoing shortwave radiation
SWRoutB						=	[];			
SWRoutB.SWRoutCeiling		=	SWRout_i(1);
SWRoutB.SWRoutWall		    =	SWRout_i(2);
SWRoutB.SWRoutInternalMass   =	SWRout_i(3);
SWRoutB.SWRoutGround		=	SWRout_i(4);


% Absorbed shortwave radiation
SWRabsB						=	[];			
SWRabsB.SWRabsCeiling		=	SWRnet_i(1);
SWRabsB.SWRabsWall		    =	SWRnet_i(2);
SWRabsB.SWRabsInternalMass   =	SWRnet_i(3);
SWRabsB.SWRabsGround        =	SWRnet_i(4);


% Energy Balance of shortwave radiation								
SWREBB				    =	[];			
SWREBB.SWREBCeiling     =	SWRinB.SWRinCeiling - SWRoutB.SWRoutCeiling - SWRabsB.SWRabsCeiling;
SWREBB.SWREBWall		=	SWRinB.SWRinWall - SWRoutB.SWRoutWall - SWRabsB.SWRabsWall;
SWREBB.SWREBInternalMass=	SWRinB.SWRinInternalMass - SWRoutB.SWRoutInternalMass - SWRabsB.SWRabsInternalMass;
SWREBB.SWREBGround      =	SWRinB.SWRinGround - SWRoutB.SWRoutGround - SWRabsB.SWRabsGround;


if abs(SWREBB.SWREBCeiling)>=10^-6 
	disp('SWREBB.SWREBCeiling is not 0. Please check SWRabsIndoors.m')
end
if abs(SWREBB.SWREBWall)>=10^-6
	disp('SWREBB.SWREBWall is not 0. Please check SWRabsIndoors.m')
end
if abs(SWREBB.SWREBInternalMass )>=10^-6
	disp('SWREBB.SWREBInternalMass	 is not 0. Please check SWRabsIndoors.m')
end
if abs(SWREBB.SWREBGround)>=10^-6
	disp('SWREBB.SWREBGround is not 0. Please check SWRabsIndoors.m')
end
