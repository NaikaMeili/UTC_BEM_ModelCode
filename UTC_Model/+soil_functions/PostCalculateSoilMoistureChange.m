function[dVdtRoofCalc,dVdtCanCalc,dVdtUrbCalc]=PostCalculateSoilMoistureChange(...
    OwaterInitial,Owater,ParSoil,FractionsRoof_Out,FractionsGround_Out,geometry_Out)

% Change in soil water volumn
%--------------------------------------------------------------------------
% Postcalculate soil water change
% Intial soil water setting in beginning of simulation (time step 1)
Vinit_Rveg   =  (OwaterInitial.OwRoofSoilVeg-ParSoil.Roof.Ohy).*ParSoil.Roof.dz;
Vinit_Gimp   =  (OwaterInitial.OwGroundSoilImp-ParSoil.Ground.Ohy).*ParSoil.Ground.dz;
Vinit_Gbare  =  (OwaterInitial.OwGroundSoilBare-ParSoil.Ground.Ohy).*ParSoil.Ground.dz;
Vinit_Gveg   =  (OwaterInitial.OwGroundSoilVeg-ParSoil.Ground.Ohy).*ParSoil.Ground.dz;

% Water volumne in soil each soil column
VRveg	=   (Owater.OwRoofSoilVeg - ParSoil.Roof.Ohy).*ParSoil.Roof.dz;
VGimp	=   (Owater.OwGroundSoilImp - ParSoil.Ground.Ohy).*ParSoil.Ground.dz;
VGbare	=   (Owater.OwGroundSoilBare - ParSoil.Ground.Ohy).*ParSoil.Ground.dz;
VGveg	=   (Owater.OwGroundSoilVeg - ParSoil.Ground.Ohy).*ParSoil.Ground.dz;

% Total water volumn in soil (canyon & roof)
VinitRoof   =   FractionsRoof_Out.fveg.*nansum(Vinit_Rveg,2);
VinitCan    =   FractionsGround_Out.fimp.*nansum(Vinit_Gimp,2)+ ...
                FractionsGround_Out.fbare.*nansum(Vinit_Gbare,2) + ...
                FractionsGround_Out.fveg.*nansum(Vinit_Gveg,2);
            
VRoof   =   FractionsRoof_Out.fveg.*sum(VRveg,2);
VCan    =   FractionsGround_Out.fimp.*sum(VGimp,2)+ ...
            FractionsGround_Out.fbare.*sum(VGbare,2) + ...
            FractionsGround_Out.fveg.*sum(VGveg,2);
        
VUrb    =   geometry_Out.wcanyon_norm.*VCan + geometry_Out.wroof_norm.*VRoof;       
        
% Change in soil volumn in soil column
VRoof_t      = VRoof;
VRoof_tm1    = [VinitRoof; VRoof(1:end-1)];
dVdtRoofCalc = VRoof_t - VRoof_tm1;

VCan_t      = VCan;
VCan_tm1    = [VinitCan; VCan(1:end-1)];
dVdtCanCalc = VCan_t - VCan_tm1;

dVdtUrbCalc     = geometry_Out.wcanyon_norm.*dVdtCanCalc + geometry_Out.wroof_norm.*dVdtRoofCalc;



