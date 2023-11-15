function[dVdtRoofCalc,dVdtCanCalc,dVdtUrbCalc]=PostCalculateSoilMoistureChange(...
    OwaterInitial,Owater,ParSoilRoof_Out,ParSoilGround_Out,FractionsRoof_Out,FractionsGround_Out,geometry_Out,ittm)

% Change in soil water volumn
%--------------------------------------------------------------------------
% Postcalculate soil water change
% Intial soil water setting in beginning of simulation (time step 1)
Vinit_Rveg   =  (OwaterInitial.OwRoofSoilVeg(:,:,ittm)-ParSoilRoof_Out(ittm).Ohy).*ParSoilRoof_Out(ittm).dz;
Vinit_Gimp   =  (OwaterInitial.OwGroundSoilImp(:,:,ittm)-ParSoilGround_Out(ittm).Ohy).*ParSoilGround_Out(ittm).dz;
Vinit_Gbare  =  (OwaterInitial.OwGroundSoilBare(:,:,ittm)-ParSoilGround_Out(ittm).Ohy).*ParSoilGround_Out(ittm).dz;
Vinit_Gveg   =  (OwaterInitial.OwGroundSoilVeg(:,:,ittm)-ParSoilGround_Out(ittm).Ohy).*ParSoilGround_Out(ittm).dz;

% Water volumne in soil each soil column
VRveg	=   (Owater.OwRoofSoilVeg(:,:,ittm) - ParSoilRoof_Out(ittm).Ohy).*ParSoilRoof_Out(ittm).dz;
VGimp	=   (Owater.OwGroundSoilImp(:,:,ittm) - ParSoilGround_Out(ittm).Ohy).*ParSoilGround_Out(ittm).dz;
VGbare	=   (Owater.OwGroundSoilBare(:,:,ittm) - ParSoilGround_Out(ittm).Ohy).*ParSoilGround_Out(ittm).dz;
VGveg	=   (Owater.OwGroundSoilVeg(:,:,ittm) - ParSoilGround_Out(ittm).Ohy).*ParSoilGround_Out(ittm).dz;

% Total water volumn in soil (canyon & roof)
VinitRoof   =   FractionsRoof_Out(ittm).fveg.*nansum(Vinit_Rveg,2);
VinitCan    =   FractionsGround_Out(ittm).fimp.*nansum(Vinit_Gimp,2)+ ...
                FractionsGround_Out(ittm).fbare.*nansum(Vinit_Gbare,2) + ...
                FractionsGround_Out(ittm).fveg.*nansum(Vinit_Gveg,2);
            
VRoof   =   FractionsRoof_Out(ittm).fveg.*sum(VRveg,2);
VCan    =   FractionsGround_Out(ittm).fimp.*sum(VGimp,2)+ ...
            FractionsGround_Out(ittm).fbare.*sum(VGbare,2) + ...
            FractionsGround_Out(ittm).fveg.*sum(VGveg,2);
        
VUrb    =   geometry_Out(ittm).wcanyon_norm.*VCan + geometry_Out(ittm).wroof_norm.*VRoof;       
        
% Change in soil volumn in soil column
VRoof_t      = VRoof;
VRoof_tm1    = [VinitRoof; VRoof(1:end-1)];
dVdtRoofCalc = VRoof_t - VRoof_tm1;

VCan_t      = VCan;
VCan_tm1    = [VinitCan; VCan(1:end-1)];
dVdtCanCalc = VCan_t - VCan_tm1;

dVdtUrbCalc     = geometry_Out(ittm).wcanyon_norm.*dVdtCanCalc + geometry_Out(ittm).wroof_norm.*dVdtRoofCalc;



