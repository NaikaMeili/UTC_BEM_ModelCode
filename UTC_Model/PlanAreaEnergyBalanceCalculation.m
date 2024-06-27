function[EnergyFluxUrban,EnergyFluxCan,EnergyFluxRoof]=PlanAreaEnergyBalanceCalculation(ViewFactor,MeteoDataRaw,...
    SWRin,SWRout,SWRabs,LWRin,LWRout,LWRabs,LEflux,Hflux,Gflux,BEMWasteHeat,...
    dStorage,GbuildInt,HbuildInt,LEbuildInt,SWRabsB,LWRabsB,...
    geometry_Out,FractionsGround_Out,PropOpticalRoof_Out,Anthropo,Gemeotry_m_Out,Figure,ittmTot,BEM_on)

for ittm=1:ittmTot
% Incoming radiation
%--------------------------------------------------------------------------
SWRin_atm			=	(MeteoDataRaw.SAB1_in + MeteoDataRaw.SAB2_in + MeteoDataRaw.SAD1_in + MeteoDataRaw.SAD2_in);
SWRinDir_atm		=	(MeteoDataRaw.SAB1_in + MeteoDataRaw.SAB2_in);
SWRinDiff_atm		=	(MeteoDataRaw.SAD1_in + MeteoDataRaw.SAD2_in);

LWRin_atm			=	MeteoDataRaw.LWR_in;

% View factors
%--------------------------------------------------------------------------
F_gs_T=ViewFactor.F_gs_T; F_gt_T=ViewFactor.F_gt_T; F_gw_T=ViewFactor.F_gw_T;
F_ww_T=ViewFactor.F_ww_T; F_wt_T=ViewFactor.F_wt_T; F_wg_T=ViewFactor.F_wg_T;
F_ws_T=ViewFactor.F_ws_T; F_sg_T=ViewFactor.F_sg_T; F_sw_T=ViewFactor.F_sw_T;
F_st_T=ViewFactor.F_st_T; F_tg_T=ViewFactor.F_tg_T; F_tw_T=ViewFactor.F_tw_T;
F_ts_T=ViewFactor.F_ts_T; F_tt_T=ViewFactor.F_tt_T;

% Normalize surface area
%--------------------------------------------------------------------------
A_s     =   geometry_Out(ittm).wcanyon; 
A_g     =   geometry_Out(ittm).wcanyon; 
A_w     =   geometry_Out(ittm).hcanyon;
A_t     =   2*2*pi*geometry_Out(ittm).radius_tree;	% There are 2 trees. Hence, the area of tree is twice a circle
fgveg   =   FractionsGround_Out(ittm).fveg; 
fgbare  =   FractionsGround_Out(ittm).fbare; 
fgimp   =   FractionsGround_Out(ittm).fimp;

% Rescaling tree absorbed radiation
SWRabs.SWRabsTree(:,1,ittm) = SWRabs.SWRabsTree(:,1,ittm).*4.*geometry_Out(ittm).radius_tree./(4.*geometry_Out(ittm).radius_tree.*pi);
LWRabs.LWRabsTree(:,1,ittm) = LWRabs.LWRabsTree(:,1,ittm).*4.*geometry_Out(ittm).radius_tree./(4.*geometry_Out(ittm).radius_tree.*pi);

% Shortwave canyon components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CanSWRin_SurfArea	=	A_g./A_g.*(SWRin.SWRinGroundVeg(:,1,ittm).*fgveg + SWRin.SWRinGroundBare(:,1,ittm).*fgbare + SWRin.SWRinGroundImp(:,1,ittm).*fgimp) + ...
						A_w./A_g.*(SWRin.SWRinWallSun(:,1,ittm) + SWRin.SWRinWallShade(:,1,ittm)) + A_t./A_g.*SWRin.SWRinTree(:,1,ittm);

CanSWRabs_SurfArea	=	A_g./A_g.*(SWRabs.SWRabsGroundVeg(:,1,ittm).*fgveg + SWRabs.SWRabsGroundBare(:,1,ittm).*fgbare + SWRabs.SWRabsGroundImp(:,1,ittm).*fgimp) + ...
						A_w./A_g.*(SWRabs.SWRabsWallSun(:,1,ittm) + SWRabs.SWRabsWallShade(:,1,ittm)) + SWRabs.SWRabsTree(:,1,ittm).*A_t./A_g;
      
CanSWRout_SurfArea	=	A_g./A_g.*(SWRout.SWRoutGroundVeg(:,1,ittm).*fgveg + SWRout.SWRoutGroundBare(:,1,ittm).*fgbare + SWRout.SWRoutGroundImp(:,1,ittm).*fgimp) +...
						A_w./A_g.*(SWRout.SWRoutWallSun(:,1,ittm) + SWRout.SWRoutWallShade(:,1,ittm)) + SWRout.SWRoutTree(:,1,ittm).*A_t./A_s;
	  
CanSWRout_Ref_to_Atm=	SWRout.SWRoutGroundVeg(:,1,ittm).*F_sg_T.*fgveg + SWRout.SWRoutGroundBare(:,1,ittm).*F_sg_T.*fgbare + SWRout.SWRoutGroundImp(:,1,ittm).*F_sg_T.*fgimp + ...
						SWRout.SWRoutWallSun(:,1,ittm).*F_sw_T + SWRout.SWRoutWallShade(:,1,ittm).*F_sw_T + SWRout.SWRoutTree(:,1,ittm).*F_st_T;
					
CanSWR_EB_SurfArea			=	CanSWRin_SurfArea - CanSWRabs_SurfArea - CanSWRout_SurfArea;
CanSWR_EB_PlanAreaCanyon	=	SWRin_atm - CanSWRabs_SurfArea - CanSWRout_Ref_to_Atm;

% Check that the sum of SWR transmitted through windows and absorbed by
% exterior surfaces is the same.
EBwallSunSWR = SWRabs.SWRabsWallSun(:,1,ittm)-SWRabs.SWRabsWallSunExt(:,1,ittm)-SWRabs.SWRabsWallSunTransmitted(:,1,ittm);
EBwallShdSWR = SWRabs.SWRabsWallShade(:,1,ittm)-SWRabs.SWRabsWallShadeExt(:,1,ittm)-SWRabs.SWRabsWallShadeTransmitted(:,1,ittm);

% Check interior SWR balance: total absorbed interior SWR should be equal
% to the transmitted radiation
SWREBinternal = Gemeotry_m_Out(ittm).Height_canyon.*(SWRabs.SWRabsWallSunTransmitted(:,1,ittm)+SWRabs.SWRabsWallShadeTransmitted(:,1,ittm))...
                -Gemeotry_m_Out(ittm).Height_canyon.*(SWRabsB.SWRabsWallsun(:,1,ittm)+SWRabsB.SWRabsWallshd(:,1,ittm)) ...
                -Gemeotry_m_Out(ittm).Width_roof.*(SWRabsB.SWRabsCeiling(:,1,ittm) + SWRabsB.SWRabsGround(:,1,ittm))...
                -Gemeotry_m_Out(ittm).Height_canyon.*SWRabsB.SWRabsInternalMass(:,1,ittm);

% Shortwave urban components
%--------------------------------------------------------------------------
% Calculate total urban including roofs
UrbanSWRin_SurfArea     =   geometry_Out(ittm).wcanyon_norm.*CanSWRin_SurfArea + geometry_Out(ittm).wroof_norm.*SWRin.SWRinTotalRoof(:,1,ittm);
UrbanSWRabs_SurfArea    =   geometry_Out(ittm).wcanyon_norm.*CanSWRabs_SurfArea + geometry_Out(ittm).wroof_norm.*SWRabs.SWRabsTotalRoof(:,1,ittm);
UrbanSWRout_SurfArea    =   geometry_Out(ittm).wcanyon_norm.*CanSWRout_SurfArea + geometry_Out(ittm).wroof_norm.*SWRout.SWRoutTotalRoof(:,1,ittm);
UrbanSWRout_Ref_to_Atm  =   geometry_Out(ittm).wcanyon_norm.*CanSWRout_Ref_to_Atm + geometry_Out(ittm).wroof_norm.*SWRout.SWRoutTotalRoof(:,1,ittm);

UrbanSWR_EB_SurfArea        =	UrbanSWRin_SurfArea - UrbanSWRabs_SurfArea - UrbanSWRout_SurfArea;
UrbanSWR_EB_PlanAreaUrban	=	SWRin_atm - UrbanSWRabs_SurfArea - UrbanSWRout_Ref_to_Atm;

% Canyon albedo calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UrbanAlbedo	=	UrbanSWRout_Ref_to_Atm./SWRin_atm;
CanAlbedo	=	CanSWRout_Ref_to_Atm./SWRin_atm;
RoofAlbedo	=	PropOpticalRoof_Out(ittm).albedo;


% Longwave canyon components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CanLWRin_SurfArea	=	LWRin.LWRinGroundVeg(:,1,ittm).*fgveg.*A_g./A_g + LWRin.LWRinGroundBare(:,1,ittm).*fgbare.*A_g./A_g + LWRin.LWRinGroundImp(:,1,ittm).*fgimp.*A_g./A_g + ...
						LWRin.LWRinWallSun(:,1,ittm).*A_w./A_g + LWRin.LWRinWallShade(:,1,ittm).*A_w./A_g + LWRin.LWRinTree(:,1,ittm).*A_t./A_g;

CanLWRabs_SurfArea	=	LWRabs.LWRabsGroundVeg(:,1,ittm).*fgveg.*A_g./A_g + LWRabs.LWRabsGroundBare(:,1,ittm).*fgbare.*A_g./A_g + LWRabs.LWRabsGroundImp(:,1,ittm)*fgimp.*A_g./A_g + ...
						LWRabs.LWRabsWallSun(:,1,ittm).*A_w./A_g + LWRabs.LWRabsWallShade(:,1,ittm).*A_w./A_g + LWRabs.LWRabsTree(:,1,ittm).*A_t./A_g;
      
CanLWRout_SurfArea	=	LWRout.LWRoutGroundVeg(:,1,ittm).*fgveg.*A_g/A_s+LWRout.LWRoutGroundBare(:,1,ittm).*fgbare.*A_g./A_s+LWRout.LWRoutGroundImp(:,1,ittm).*fgimp.*A_g./A_s+...
						LWRout.LWRoutWallSun(:,1,ittm).*A_w./A_s+LWRout.LWRoutWallShade(:,1,ittm).*A_w./A_s+LWRout.LWRoutTree(:,1,ittm).*A_t./A_s;
	  
CanLWRout_Ref_to_Atm=	LWRout.LWRoutGroundVeg(:,1,ittm).*F_sg_T.*fgveg + LWRout.LWRoutGroundBare(:,1,ittm).*F_sg_T.*fgbare + LWRout.LWRoutGroundImp(:,1,ittm).*F_sg_T.*fgimp + ...
						LWRout.LWRoutWallSun(:,1,ittm).*F_sw_T + LWRout.LWRoutWallShade(:,1,ittm).*F_sw_T + LWRout.LWRoutTree(:,1,ittm).*F_st_T;
					
CanLWR_EB_SurfArea			=	CanLWRin_SurfArea - CanLWRabs_SurfArea - CanLWRout_SurfArea;
CanLWR_EB_PlanAreaCanyon	=	LWRin_atm - CanLWRabs_SurfArea - CanLWRout_Ref_to_Atm;


% Longwave urban components
%--------------------------------------------------------------------------
% Calculate total urban including roofs
UrbanLWRin_SurfArea     =   geometry_Out(ittm).wcanyon_norm.*CanLWRin_SurfArea + geometry_Out(ittm).wroof_norm.*LWRin.LWRinTotalRoof(:,1,ittm);
UrbanLWRabs_SurfArea    =   geometry_Out(ittm).wcanyon_norm.*CanLWRabs_SurfArea + geometry_Out(ittm).wroof_norm.*LWRabs.LWRabsTotalRoof(:,1,ittm);
UrbanLWRout_SurfArea    =   geometry_Out(ittm).wcanyon_norm.*CanLWRout_SurfArea + geometry_Out(ittm).wroof_norm.*LWRout.LWRoutTotalRoof(:,1,ittm);
UrbanLWRout_Ref_to_Atm  =   geometry_Out(ittm).wcanyon_norm.*CanLWRout_Ref_to_Atm + geometry_Out(ittm).wroof_norm.*LWRout.LWRoutTotalRoof(:,1,ittm);

UrbanLWR_EB_SurfArea        =	UrbanLWRin_SurfArea - UrbanLWRabs_SurfArea - UrbanLWRout_SurfArea;
UrbanLWR_EB_PlanAreaUrban	=	LWRin_atm - UrbanLWRabs_SurfArea - UrbanLWRout_Ref_to_Atm;


% Urban energy balance components
%--------------------------------------------------------------------------
SWRabs_Urban    = SWRabs.SWRabsTotalUrban(:,1,ittm);
LWRabs_Urban    = LWRabs.LWRabsTotalUrban(:,1,ittm);
LE_Urban        = LEflux.LEfluxUrban(:,1,ittm);
H_Urban         = Hflux.HfluxUrban(:,1,ittm);
if BEM_on==1
    Gground_Urban   = geometry_Out(ittm).wcanyon_norm.*Gflux.G1Ground(:,1,ittm) + geometry_Out(ittm).wroof_norm.*GbuildInt.Gfloor(:,1,ittm);
    % Change in heat storage in the building envelope
    dSdt_buildEnv   = geometry_Out(ittm).wcanyon_norm.*(A_w./A_g.*(dStorage.dsWallSun(:,1,ittm)+dStorage.dsWallShade(:,1,ittm)))...
                      + geometry_Out(ittm).wroof_norm.*(dStorage.dsRoof(:,1,ittm))...
                      + Gemeotry_m_Out(ittm).Height_canyon./(Gemeotry_m_Out(ittm).Width_canyon + Gemeotry_m_Out(ittm).Width_roof).*GbuildInt.dSinternalMass(:,1,ittm);
    G1Building = 0;
else
    Gground_Urban   = geometry_Out(ittm).wcanyon_norm.*Gflux.G1Ground(:,1,ittm);
    % Change in heat storage in the building envelope
    G1Building   = geometry_Out(ittm).wcanyon_norm.*(A_w./A_g.*(Gflux.G1WallSun(:,1,ittm)+Gflux.G1WallShade(:,1,ittm)))...
                      + geometry_Out(ittm).wroof_norm.*(Gflux.G1Roof(:,1,ittm));
    dSdt_buildEnv = 0;
end
% Change in heat storage in the air
dSdt_Air        = geometry_Out(ittm).wcanyon_norm.*(Hflux.dS_H_air(:,1,ittm) + LEflux.dS_LE_air(:,1,ittm))...
                  +geometry_Out(ittm).wroof_norm.*(HbuildInt.dSH_air(:,1,ittm) + LEbuildInt.dSLE_air(:,1,ittm));

% Total anthropogenic heat input
Qanth           = Anthropo.Qf_canyon(:,1,ittm).*geometry_Out(ittm).wcanyon_norm + Anthropo.Qf_roof(:,1,ittm).*geometry_Out(ittm).wroof_norm...
                    + BEMWasteHeat.TotAnthInput_URB(:,1,ittm);

% Anthropogenic heat input due to condensation of water during AC process
QanthACcondensation =  BEMWasteHeat.WaterFromAC_Can.*geometry_Out(ittm).wcanyon_norm;

% Total urban energy balance
if BEM_on==1
    EBtot = SWRabs_Urban + LWRabs_Urban + Qanth + QanthACcondensation...
            - LE_Urban - H_Urban - Gground_Urban - dSdt_buildEnv - dSdt_Air;
else
    EBtot = SWRabs_Urban + LWRabs_Urban + Qanth...
            - LE_Urban - H_Urban - Gground_Urban - G1Building - dSdt_Air;
end


% figure
% plot(MeteoDataRaw.Date,EBtot); title(['EB urban, ittm = ' num2str(ittm)])

% Critical Energy Balance components out
%--------------------------------------------------------------------------
% Urban
%--------------------------------------------------------------------------
EnergyFluxUrban.SWRin_PlanArea(:,1,ittm)          = SWRin_atm;
EnergyFluxUrban.SWRin_SurfArea(:,1,ittm)          = UrbanSWRin_SurfArea;
EnergyFluxUrban.SWRabs_SurfArea(:,1,ittm)         = UrbanSWRabs_SurfArea;
EnergyFluxUrban.SWRout_SurfArea(:,1,ittm)         = UrbanSWRout_SurfArea;
EnergyFluxUrban.SWRout_Ref_to_Atm(:,1,ittm)       = UrbanSWRout_Ref_to_Atm;
EnergyFluxUrban.SWREB_SurfArea(:,1,ittm)          = UrbanSWR_EB_SurfArea;
EnergyFluxUrban.SWREB_PlanAreaUrban(:,1,ittm)     = UrbanSWR_EB_PlanAreaUrban;

EnergyFluxUrban.LWRin_PlanArea(:,1,ittm)          = LWRin_atm;
EnergyFluxUrban.LWRin_SurfArea(:,1,ittm)         = UrbanLWRin_SurfArea;
EnergyFluxUrban.LWRabs_SurfArea(:,1,ittm)         = UrbanLWRabs_SurfArea;
EnergyFluxUrban.LWRout_SurfArea(:,1,ittm)         = UrbanLWRout_SurfArea;
EnergyFluxUrban.LWRout_Ref_to_Atm(:,1,ittm)       = UrbanLWRout_Ref_to_Atm;
EnergyFluxUrban.LWREB_SurfArea(:,1,ittm)          = UrbanLWR_EB_SurfArea;
EnergyFluxUrban.LWREB_PlanAreaUrban(:,1,ittm)     = UrbanLWR_EB_PlanAreaUrban;

EnergyFluxUrban.UrbanAlbedo(:,1,ittm)  = UrbanAlbedo;

EnergyFluxUrban.SWRabs(:,1,ittm)	        = SWRabs_Urban;
EnergyFluxUrban.LWRabs(:,1,ittm)	        = LWRabs_Urban;
EnergyFluxUrban.LEflux(:,1,ittm)	        = LE_Urban;
EnergyFluxUrban.Hflux(:,1,ittm)	            = H_Urban;
EnergyFluxUrban.GfluxGround(:,1,ittm)	    = Gground_Urban;
EnergyFluxUrban.G1Building(:,1,ittm)        = G1Building;
EnergyFluxUrban.dSdtBuildEnv(:,1,ittm)      = dSdt_buildEnv;
EnergyFluxUrban.dSdt_Air(:,1,ittm)          = dSdt_Air;
EnergyFluxUrban.Qanth(:,1,ittm)             = Qanth;
EnergyFluxUrban.QanthACcondensation(:,1,ittm)= QanthACcondensation;
EnergyFluxUrban.EB(:,1,ittm)                = EBtot;

% Canyon
%--------------------------------------------------------------------------
EnergyFluxCan.SWRin_SurfArea(:,1,ittm)          = CanSWRin_SurfArea;
EnergyFluxCan.SWRabs_SurfArea(:,1,ittm)         = CanSWRabs_SurfArea;
EnergyFluxCan.SWRout_SurfArea(:,1,ittm)         = CanSWRout_SurfArea;
EnergyFluxCan.SWRout_Ref_to_Atm(:,1,ittm)       = CanSWRout_Ref_to_Atm;
EnergyFluxCan.SWREB_SurfArea(:,1,ittm)          = CanSWR_EB_SurfArea;
EnergyFluxCan.SWREB_PlanAreaCanyon(:,1,ittm)	= CanSWR_EB_PlanAreaCanyon;

EnergyFluxCan.LWRin_SurfArea(:,1,ittm)          = CanLWRin_SurfArea;
EnergyFluxCan.LWRabs_SurfArea(:,1,ittm)         = CanLWRabs_SurfArea;
EnergyFluxCan.LWRout_SurfArea(:,1,ittm)         = CanLWRout_SurfArea;
EnergyFluxCan.LWRout_Ref_to_Atm(:,1,ittm)       = CanLWRout_Ref_to_Atm;
EnergyFluxCan.LWREB_SurfArea(:,1,ittm)          = CanLWR_EB_SurfArea;
EnergyFluxCan.LWREB_PlanAreaCanyon(:,1,ittm)  = CanLWR_EB_PlanAreaCanyon;

EnergyFluxCan.CanAlbedo(:,1,ittm)  = CanAlbedo;

EnergyFluxCan.SWRabs(:,1,ittm)	= SWRabs.SWRabsTotalCanyon(:,1,ittm)-A_w./A_g.*(SWRabs.SWRabsWallSunTransmitted(:,1,ittm)+SWRabs.SWRabsWallShadeTransmitted(:,1,ittm));
EnergyFluxCan.LWRabs(:,1,ittm)	= LWRabs.LWRabsTotalCanyon(:,1,ittm);
EnergyFluxCan.LEflux(:,1,ittm)	= LEflux.LEfluxCanyon(:,1,ittm);
EnergyFluxCan.Hflux(:,1,ittm)     = Hflux.HfluxCanyon(:,1,ittm);
EnergyFluxCan.Gflux(:,1,ittm)     = Gflux.G1Canyon(:,1,ittm);
EnergyFluxCan.dSdt_Air(:,1,ittm)  = (Hflux.dS_H_air(:,1,ittm) + LEflux.dS_LE_air(:,1,ittm));
EnergyFluxCan.Qanth(:,1,ittm)     = Anthropo.Qf_canyon(:,1,ittm) + ...
                                        + BEMWasteHeat.SensibleFromVent_Can(:,1,ittm) + BEMWasteHeat.SensibleFromAC_Can(:,1,ittm) + BEMWasteHeat.SensibleFromHeat_Can(:,1,ittm)...
                                        + BEMWasteHeat.LatentFromVent_Can(:,1,ittm) + BEMWasteHeat.LatentFromAC_Can(:,1,ittm) +  BEMWasteHeat.LatentFromHeat_Can(:,1,ittm);

% Anthropogenic heat input due to condensation of water during AC process
EnergyFluxCan.QanthACcondensation(:,1,ittm) =  BEMWasteHeat.WaterFromAC_Can(:,1,ittm);

EnergyFluxCan.EB(:,1,ittm)        = EnergyFluxCan.SWRabs(:,1,ittm) + EnergyFluxCan.LWRabs(:,1,ittm) + EnergyFluxCan.Qanth(:,1,ittm) + EnergyFluxCan.QanthACcondensation(:,1,ittm)...
                            - EnergyFluxCan.LEflux(:,1,ittm) - EnergyFluxCan.Hflux(:,1,ittm) - EnergyFluxCan.Gflux(:,1,ittm) - EnergyFluxCan.dSdt_Air(:,1,ittm);


% figure
% plot(EnergyFluxCan.EB(:,1,ittm)); title(['EB canyon, ittm = ' num2str(ittm)])


% Roof
%--------------------------------------------------------------------------
EnergyFluxRoof.SWRin(:,1,ittm)	= SWRin.SWRinTotalRoof(:,1,ittm);
EnergyFluxRoof.SWRout(:,1,ittm)	= SWRout.SWRoutTotalRoof(:,1,ittm);
EnergyFluxRoof.LWRin(:,1,ittm)	= LWRin.LWRinTotalRoof(:,1,ittm);
EnergyFluxRoof.LWRout(:,1,ittm)	= LWRout.LWRoutTotalRoof(:,1,ittm);

EnergyFluxRoof.RoofAlbedo(:,1,ittm)  = RoofAlbedo;

EnergyFluxRoof.SWRabs(:,1,ittm)	= SWRabs.SWRabsTotalRoof(:,1,ittm);
EnergyFluxRoof.LWRabs(:,1,ittm)	= LWRabs.LWRabsTotalRoof(:,1,ittm);
EnergyFluxRoof.LEflux(:,1,ittm)	= LEflux.LEfluxRoof(:,1,ittm);
EnergyFluxRoof.Hflux(:,1,ittm)   = Hflux.HfluxRoof(:,1,ittm);
EnergyFluxRoof.Gflux(:,1,ittm)    = Gflux.G1Roof(:,1,ittm);
EnergyFluxRoof.Qanth(:,1,ittm)    = Anthropo.Qf_roof(:,1,ittm);
EnergyFluxRoof.EB(:,1,ittm)       = EnergyFluxRoof.SWRabs(:,1,ittm) + EnergyFluxRoof.LWRabs(:,1,ittm) - EnergyFluxRoof.LEflux(:,1,ittm) - ...
                            EnergyFluxRoof.Hflux(:,1,ittm) - EnergyFluxRoof.Gflux(:,1,ittm) + EnergyFluxRoof.Qanth(:,1,ittm);

EnergyFluxRoof.SWREB(:,1,ittm)	= SWRin.SWRinTotalRoof(:,1,ittm) - SWRabs.SWRabsTotalRoof(:,1,ittm) - SWRout.SWRoutTotalRoof(:,1,ittm);
EnergyFluxRoof.LWREB(:,1,ittm)    = LWRin.LWRinTotalRoof(:,1,ittm) - LWRabs.LWRabsTotalRoof(:,1,ittm) - LWRout.LWRoutTotalRoof(:,1,ittm);

% figure
% plot(EnergyFluxRoof.EB(:,1,ittm)) ; title('EB roof')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot graphs
if Figure==1
for ittm=1:ittmTot
    
TTUrban = table(hour(MeteoDataRaw.Date),month(MeteoDataRaw.Date),...
            EnergyFluxUrban.SWRin_PlanArea(:,1,ittm),EnergyFluxUrban.SWRabs_SurfArea(:,1,ittm),EnergyFluxUrban.SWRout_Ref_to_Atm(:,1,ittm),...
            EnergyFluxUrban.LWRin_PlanArea(:,1,ittm),EnergyFluxUrban.LWRabs_SurfArea(:,1,ittm),EnergyFluxUrban.LWRout_Ref_to_Atm(:,1,ittm),...
            EnergyFluxUrban.LEflux(:,1,ittm),EnergyFluxUrban.Hflux(:,1,ittm),EnergyFluxUrban.GfluxGround(:,1,ittm) + EnergyFluxUrban.dSdtBuildEnv(:,1,ittm) + EnergyFluxUrban.G1Building(:,1,ittm),...
            EnergyFluxUrban.Qanth(:,1,ittm),EnergyFluxUrban.QanthACcondensation(:,1,ittm),EnergyFluxUrban.EB(:,1,ittm),EnergyFluxUrban.UrbanAlbedo(:,1,ittm),...
            EnergyFluxUrban.Hflux(:,1,ittm)./EnergyFluxUrban.LEflux(:,1,ittm));

TTUrban.Properties.VariableNames = {'Hour','Month','SWRin','SWRabs','SWRout','LWRin','LWRabs','LWRout',...
    'LE','H','G','Qanth','QanthACcondensation','EB','Albedo','BowenRatio'};

TTUrbanDiurnal = varfun(@nanmean,TTUrban,'GroupingVariables','Hour');
TTUrbanSeasonal = varfun(@nanmean,TTUrban,'GroupingVariables','Month');

TTUrbanDiurnalMedian = varfun(@median,TTUrban,'GroupingVariables','Hour');
TTUrbanSeasonalMedian = varfun(@median,TTUrban,'GroupingVariables','Month');


% Energy budget
f1 = figure;
set(f1, 'Units','centimeters','Position', [1 1 19 7])

t = tiledlayout(1,3);
t.Padding = 'compact'; %t.TileSpacing = 'compact';

nexttile
plot(MeteoDataRaw.Date,EnergyFluxUrban.EB(:,1,ittm),'k','DisplayName','EB')
xlabel('time'); ylabel('EB W/m^{2}'); title('Time series');
%legend

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_EB,'k','LineWidth',1.5,'DisplayName','EB')
xlim([0 23]); xlabel('hour'); ylabel('EB W/m^{2}'); title('Diurnal');
%legend

nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_EB,'k','LineWidth',1.5,'DisplayName','EB')
xlim([1 12]); xlabel('Month'); ylabel('EB W/m^{2}'); title('Seasonal');
%legend
sgtitle(['Energy budget closure ittm = ' num2str(ittm),', EB is not fully closed if parital AC is applied'])

% Albedo und Bowen ratio
f1 = figure;
set(f1, 'Units','centimeters','Position', [1 1 14 10])

t = tiledlayout(2,2);
t.Padding = 'compact'; %t.TileSpacing = 'compact';

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_Albedo,'k','LineWidth',1.5,'DisplayName','albedo')
xlim([0 23]); xlabel('hour'); ylabel('Albedo (-)'); title('Albedo = SWR_{out}/SWR_{in}');
subtitle('Diurnal');
%legend

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnalMedian.median_BowenRatio,'k','LineWidth',1.5,'DisplayName','Bowen ratio')
xlim([0 23]); xlabel('hour'); ylabel('Bowen Ratio (-)'); title('Bowen ratio = H/LE, median');
subtitle('Diurnal');
%legend

nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_Albedo,'k','LineWidth',1.5,'DisplayName','albedo')
xlim([1 12]); xlabel('month'); ylabel('Albedo (-)'); subtitle('Seasonal');
%legend

nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonalMedian.median_BowenRatio,'k','LineWidth',1.5,'DisplayName','Bowen ratio')
xlim([1 12]); xlabel('month'); ylabel('Bowen Ratio (-)'); subtitle('Seasonal');
%legend

% Energy Fluxes
f1 = figure;
set(f1, 'Units','centimeters','Position', [1 1 19 14])

t = tiledlayout(2,3);
t.Padding = 'compact'; %t.TileSpacing = 'compact';

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_SWRin,'k','LineWidth',1.5,'DisplayName','SWR_{in}')
hold on
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_SWRabs,'r','LineWidth',1.5,'DisplayName','SWR_{abs}')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_SWRout,'b','LineWidth',1.5,'DisplayName','SWR_{out}')
xlim([0 23]); xlabel('hour'); ylabel('SWR W/m^{2}'); title('SWR'); subtitle('Diurnal');
legend('Location','southoutside','NumColumns',2)

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_LWRin,'k','LineWidth',1.5,'DisplayName','LWR_{in}')
hold on
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_LWRabs,'r','LineWidth',1.5,'DisplayName','LWR_{abs}')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_LWRout,'b','LineWidth',1.5,'DisplayName','LWR_{out}')
xlim([0 23]); xlabel('hour'); ylabel('LWR W/m^{2}'); title('LWR'); subtitle('Diurnal');
legend('Location','southoutside','NumColumns',2)

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_SWRabs,'k','LineWidth',1.5,'DisplayName','SWR_{abs}')
hold on
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_H,'r','LineWidth',1.5,'DisplayName','H')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_LE,'g','LineWidth',1.5,'DisplayName','LE')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_G,'b','LineWidth',1.5,'DisplayName','G')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_LWRabs,'m','LineWidth',1.5,'DisplayName','LWR_{abs}')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_Qanth,'k:','LineWidth',1.5,'DisplayName','Q_{f}')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_QanthACcondensation,'k--','LineWidth',1.5,'DisplayName','Q_{f,AC,cond}')
xlim([0 23]); xlabel('hour'); ylabel('EB W/m^{2}'); title('Energy fluxes'); subtitle('Diurnal');
legend('Location','southoutside','NumColumns',2)

nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_SWRin,'k','LineWidth',1.5,'DisplayName','SWR_{in}')
hold on
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_SWRabs,'r','LineWidth',1.5,'DisplayName','SWR_{abs}')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_SWRout,'b','LineWidth',1.5,'DisplayName','SWR_{out}')
xlim([1 12]); xlabel('month'); ylabel('SWR W/m^{2}'); subtitle('Seasonal'); %title('SWR'); 
%legend

nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_LWRin,'k','LineWidth',1.5,'DisplayName','LWR_{in}')
hold on
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_LWRabs,'r','LineWidth',1.5,'DisplayName','LWR_{abs}')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_LWRout,'b','LineWidth',1.5,'DisplayName','LWR_{out}')
xlim([1 12]); xlabel('month'); ylabel('LWR W/m^{2}');  subtitle('Seasonal'); %title('LWR');
%legend

nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_SWRabs,'k','LineWidth',1.5,'DisplayName','SWR_{abs}')
hold on
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_H,'r','LineWidth',1.5,'DisplayName','H')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_LE,'g','LineWidth',1.5,'DisplayName','LE')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_G,'b','LineWidth',1.5,'DisplayName','G')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_LWRabs,'m','LineWidth',1.5,'DisplayName','LWR_{abs}')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_Qanth,'k:','LineWidth',1.5,'DisplayName','Q_{f}')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_QanthACcondensation,'k--','LineWidth',1.5,'DisplayName','Q_{f,AC,cond}')
xlim([1 12]); xlabel('month'); ylabel('EB W/m^{2}'); subtitle('Seasonal'); %title('EB');
%legend
sgtitle(['Energy fluxes, ittm = ' num2str(ittm)])

end
end
