function[LWRinB,LWRoutB,LWRabsB]=LWRabsIndoors(Tinwallsun,Tinwallshd,Tceiling,Tground,Tintmass,Hbuild,Wroof,ec,eg,em,ew)


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


[LWRinB_wsun,LWRoutB_wsun,LWRabsB_wsun]=BuildingEnergyModel.LWRabsBuildingHalf(Tceiling,Tinwallsun,Tintmass,Tground,...
    A_c,A_g,A_h,F_gc,F_gw,F_ww,F_wg,F_wc,F_cg,F_cw,ec,eg,ew,em);

[LWRinB_wshd,LWRoutB_wshd,LWRabsB_wshd]=BuildingEnergyModel.LWRabsBuildingHalf(Tceiling,Tinwallshd,Tintmass,Tground,...
    A_c,A_g,A_h,F_gc,F_gw,F_ww,F_wg,F_wc,F_cg,F_cw,ec,eg,ew,em);


% Longwave radiation by each surface per m^2 surface area
% Incoming longwave radiation		
LWRinB.LWRinCeiling     = (LWRinB_wsun.LWRinCeiling + LWRinB_wshd.LWRinCeiling)./2;
LWRinB.LWRinWallsun     = LWRinB_wsun.LWRinWall;
LWRinB.LWRinWallshd     = LWRinB_wshd.LWRinWall;
LWRinB.LWRinInternalMass= LWRinB_wsun.LWRinInternalMass + LWRinB_wsun.LWRinInternalMass;
LWRinB.LWRinGround      = (LWRinB_wsun.LWRinGround + LWRinB_wshd.LWRinGround)./2;
								
% Outgoing longwave radiation	
LWRoutB.LWRoutCeiling       = (LWRoutB_wsun.LWRoutCeiling + LWRoutB_wshd.LWRoutCeiling)./2;
LWRoutB.LWRoutWallsun       = LWRoutB_wsun.LWRoutWall;
LWRoutB.LWRoutWallshd       = LWRoutB_wshd.LWRoutWall;
LWRoutB.LWRoutInternalMass  = LWRoutB_wsun.LWRoutInternalMass + LWRoutB_wshd.LWRoutInternalMass;
LWRoutB.LWRoutGround		= (LWRoutB_wsun.LWRoutGround + LWRoutB_wshd.LWRoutGround)./2;

% Absorbed longwave radiation	
LWRabsB.LWRabsCeiling		= (LWRabsB_wsun.LWRabsCeiling + LWRabsB_wshd.LWRabsCeiling)./2;
LWRabsB.LWRabsWallsun       = LWRabsB_wsun.LWRabsWall;
LWRabsB.LWRabsWallshd       = LWRabsB_wshd.LWRabsWall;
LWRabsB.LWRabsInternalMass  = LWRabsB_wsun.LWRabsInternalMass + LWRabsB_wshd.LWRabsInternalMass;
LWRabsB.LWRabsGround        = (LWRabsB_wsun.LWRabsGround + LWRabsB_wshd.LWRabsGround)./2;



