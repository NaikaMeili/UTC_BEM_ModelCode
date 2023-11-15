function[SWRinB,SWRoutB,SWRabsB]=SWRabsIndoors(SWRinWsun,SWRinWshd,Hbuild,Wroof,abc,abw,abg,abm)


% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Due to internal wall in the middle of the building, the roof area is
% halfed for the calculation of diffuse radiation reflection
Wroofint = Wroof./2;

% Calcluate view factors in the building interior
[F_gc,F_gw,F_ww,F_wg,F_wc,F_cg,F_cw,~]=BuildingEnergyModel.ViewFactorInternal(Hbuild,Wroofint);


% normalized surface areas
A_c		=	Wroofint./Wroofint;
A_g		=	Wroofint./Wroofint;
A_h		=	Hbuild./Wroofint;

% Check if view factors add up to 1
SVF(1) = F_gc + 2*F_gw; SVF(2) = F_ww + F_wg + F_wc; SVF(3) = F_cg + 2*F_cw;
SVF2(1) = F_gc + 2*F_wc*A_h; SVF2(2) = F_cg + 2*F_wg*A_h; SVF2(3) = F_ww + F_cw/A_h + F_gw/A_h;

for i=length(SVF)
	if SVF(i)<0.999 || SVF(i)>1.001
	disp('The view factors do not add up to 1 for a canyon with trees')
	end
end



% Calculate SWR absorbed in both parts of the building
[SWRinB_wsun,SWRoutB_wsun,SWRabsB_wsun]=BuildingEnergyModel.SWRabsBuildingHalf(A_c,A_g,A_h,SWRinWsun,...
    F_gc,F_gw,F_ww,F_wg,F_wc,F_cg,F_cw,abc,abw,abm,abg);

[SWRinB_wshd,SWRoutB_wshd,SWRabsB_wshd]=BuildingEnergyModel.SWRabsBuildingHalf(A_c,A_g,A_h,SWRinWshd,...
    F_gc,F_gw,F_ww,F_wg,F_wc,F_cg,F_cw,abc,abw,abm,abg);


% Shortwave radiation by each surface per m^2 surface area
% Incoming shortwave radiation		
SWRinB.SWRinCeiling     = (SWRinB_wsun.SWRinCeiling + SWRinB_wshd.SWRinCeiling)./2;
SWRinB.SWRinWallsun     = SWRinB_wsun.SWRinWall;
SWRinB.SWRinWallshd     = SWRinB_wshd.SWRinWall;
SWRinB.SWRinInternalMass= SWRinB_wsun.SWRinInternalMass + SWRinB_wsun.SWRinInternalMass;
SWRinB.SWRinGround      = (SWRinB_wsun.SWRinGround + SWRinB_wshd.SWRinGround)./2;
								
% Outgoing shortwave radiation	
SWRoutB.SWRoutCeiling       = (SWRoutB_wsun.SWRoutCeiling + SWRoutB_wshd.SWRoutCeiling)./2;
SWRoutB.SWRoutWallsun       = SWRoutB_wsun.SWRoutWall;
SWRoutB.SWRoutWallshd       = SWRoutB_wshd.SWRoutWall;
SWRoutB.SWRoutInternalMass  = SWRoutB_wsun.SWRoutInternalMass + SWRoutB_wshd.SWRoutInternalMass;
SWRoutB.SWRoutGround		= (SWRoutB_wsun.SWRoutGround + SWRoutB_wshd.SWRoutGround)./2;

% Absorbed shortwave radiation	
SWRabsB.SWRabsCeiling		= (SWRabsB_wsun.SWRabsCeiling + SWRabsB_wshd.SWRabsCeiling)./2;
SWRabsB.SWRabsWallsun       = SWRabsB_wsun.SWRabsWall;
SWRabsB.SWRabsWallshd       = SWRabsB_wshd.SWRabsWall;
SWRabsB.SWRabsInternalMass  = SWRabsB_wsun.SWRabsInternalMass + SWRabsB_wshd.SWRabsInternalMass;
SWRabsB.SWRabsGround        = (SWRabsB_wsun.SWRabsGround + SWRabsB_wshd.SWRabsGround)./2;


% Energy blanace check
SWREBinternal = Hbuild.*(SWRinWsun+SWRinWshd) - Wroof.*SWRabsB.SWRabsCeiling...
    - Hbuild.*(SWRabsB.SWRabsWallsun+SWRabsB.SWRabsWallshd+SWRabsB.SWRabsInternalMass)...
    - Wroof.*SWRabsB.SWRabsGround;


if abs(SWREBinternal)>=10^-6 
	disp('Building interior shortwave radiation balance is not 0. Please check SWRabsIndoors.m')
end


