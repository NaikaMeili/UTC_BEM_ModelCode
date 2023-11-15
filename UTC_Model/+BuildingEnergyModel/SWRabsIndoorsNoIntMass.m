function[SWRinB,SWRoutB,SWRabsB,SWREBB]=SWRabsIndoorsNoIntMass(SWRinWsun,SWRinWshd,Hbuild,Wroof,abc,abw,abg)


% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Due to internal wall in the middle of the building, the roof area is
% halfed for the calculation of diffuse radiation reflection
% Wroof = Wroof;

% Calcluate view factors in the building interior
[F_gc,F_gw,F_ww,F_wg,F_wc,F_cg,F_cw,~]=BuildingEnergyModel.ViewFactorInternal(Hbuild,Wroof);


% normalized surface areas
A_c		=	Wroof./Wroof;
A_g		=	Wroof./Wroof;
A_h		=	Hbuild./Wroof;

% Check if view factors add up to 1
SVF(1) = F_gc + 2*F_gw; SVF(2) = F_ww + F_wg + F_wc; SVF(3) = F_cg + 2*F_cw;
SVF2(1) = F_gc + 2*F_wc*A_h; SVF2(2) = F_cg + 2*F_wg*A_h; SVF2(3) = F_ww + F_cw/A_h + F_gw/A_h;

for i=length(SVF)
	if SVF(i)<0.999 || SVF(i)>1.001
	disp('The view factors do not add up to 1 for a canyon with trees')
	end
end


% Calculations for infinite reflections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sequence in Vectors :
% Ceiling
% Sunlit wall
% Shaded wall
% Ground
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Albedos
ai	=	[abc;abw;abw;abg];

% View factor matrix to solve for infinite reflections equation
% Omega_i = Tij*Bj
% B_i		=	[Bceiling; Bwallsun; Bwallshade; Bground];
Tij	=	[1, -abc*F_cw, -abc*F_cw, -abc*F_cg;...
		 -abw*F_wc, 1, -abw*F_ww, -abw*F_wg;...
         -abw*F_wc, -abw*F_ww, 1, -abw*F_wc;...
         -abg*F_gc, -abg*F_gw, -abg*F_gw, 1];
	
% Incoming shortwave radiation from sky
Omega_i	=	[abc*(F_cw*SWRinWsun + F_cw*SWRinWshd);...
			abw*F_ww*SWRinWshd;...
			abw*F_ww*SWRinWsun;...
			abg*(F_gw*SWRinWsun + F_gw*SWRinWshd)];
	
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
		
SWRdir_i=	[F_cw*SWRinWsun + F_cw*SWRinWshd;...
			F_ww*SWRinWshd;...
			F_ww*SWRinWsun;...
			F_gw*SWRinWsun + F_gw*SWRinWshd];

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

EBindoors = A_h/A_g.*(SWRinWsun+SWRinWshd) - TotalSWRSurface_abs;

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
SWRinB.SWRinWallsun     =	SWRin_i(2);
SWRinB.SWRinWallshd   =	SWRin_i(3);
SWRinB.SWRinGround      =	SWRin_i(4);

								
% Outgoing shortwave radiation
SWRoutB						=	[];			
SWRoutB.SWRoutCeiling		=	SWRout_i(1);
SWRoutB.SWRoutWallsun		=	SWRout_i(2);
SWRoutB.SWRoutWallshd		=	SWRout_i(3);
SWRoutB.SWRoutGround		=	SWRout_i(4);


% Absorbed shortwave radiation
SWRabsB						=	[];			
SWRabsB.SWRabsCeiling		=	SWRnet_i(1);
SWRabsB.SWRabsWallsun		=	SWRnet_i(2);
SWRabsB.SWRabsWallshd		=	SWRnet_i(3);
SWRabsB.SWRabsGround       =	SWRnet_i(4);


% Energy Balance of shortwave radiation								
SWREBB				    =	[];			
SWREBB.SWREBCeiling     =	SWRinB.SWRinCeiling - SWRoutB.SWRoutCeiling - SWRabsB.SWRabsCeiling;
SWREBB.SWREBWallsun		=	SWRinB.SWRinWallsun - SWRoutB.SWRoutWallsun - SWRabsB.SWRabsWallsun;
SWREBB.SWREBWallshd   =	SWRinB.SWRinWallshd - SWRoutB.SWRoutWallshd - SWRabsB.SWRabsWallshd;
SWREBB.SWREBGround      =	SWRinB.SWRinGround - SWRoutB.SWRoutGround - SWRabsB.SWRabsGround;


if abs(SWREBB.SWREBCeiling)>=10^-6 
	disp('SWREBB.SWREBCeiling is not 0. Please check SWRabsIndoors.m')
end
if abs(SWREBB.SWREBWallsun)>=10^-6
	disp('SWREBB.SWREBWallsun is not 0. Please check SWRabsIndoors.m')
end
if abs(SWREBB.SWREBWallshd )>=10^-6
	disp('SWREBB.SWREBWallshade	 is not 0. Please check SWRabsIndoors.m')
end
if abs(SWREBB.SWREBGround)>=10^-6
	disp('SWREBB.SWREBGround is not 0. Please check SWRabsIndoors.m')
end



