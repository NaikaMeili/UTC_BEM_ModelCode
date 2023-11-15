function[fconv,rsRoofPreCalc,rsGroundPreCalc,rsTreePreCalc]=PrecalculateForFasterNumericalSolution(ittn,ittm,...
         TempVec_ittm,Humidity_ittm,ParVegGround,SoilPotW_ittm,CiCO2Leaf_ittm,...
         MeteoData,HumidityAtm,geometry,FractionsGround,ParTree,...
		 PropOpticalGround,PropOpticalWall,PropOpticalTree,ParVegTree,...
         SunPosition,ViewFactor,ParWindows,BEM_on,ParVegRoof,...
         PropOpticalRoof,FractionsRoof,RES,Gemeotry_m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate enhancement factor based on ra_original of previous time step
if ittn==1
    fconv = 0;
else
    if TempVec_ittm.TCanyon-TempVec_ittm.Tatm > 0.1
    hPBL = 1000;

    [dcan,zomcan,~,~,~,~]=resistance_functions.WindProfile_Canyon(...
	Gemeotry_m.Height_canyon,Gemeotry_m.Height_tree,Gemeotry_m.Radius_tree,...
    Gemeotry_m.Width_canyon,Gemeotry_m.Width_roof,ParVegTree.Kopt,ParVegTree.LAI,...
    MeteoData.Zatm,MeteoData.Uatm,Gemeotry_m.Height_canyon,ParTree.trees,1.5,0.01,...
    Gemeotry_m.Hcan_max,Gemeotry_m.Hcan_std);

    zom_town	=	zomcan;			% Momentum roughness length of canyon, calculated according to McDonald 1998
    zoh_town	=	zom_town/10;	% Heat roughness length of canyon

    [fconv,~,~,~]=resistance_functions.EnhancementFactorRaPleim(RES.raCanyontoAtmOrig(ittn-1,1,ittm),zom_town,zoh_town,dcan,MeteoData.Zatm,MeteoData.Uatm,hPBL);
    
    else
        fconv = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Precalculate stomatal resistance for roof
if (ParVegRoof.LAI > 0 && FractionsRoof.fveg>0) 
    if ittn==1
        ra = 100; rb = 50;
        % Stomatal resistance,roughly 200-300 during the day and ca 3000 during the night    
        [rs_sun,rs_shd,Ci_sun,Ci_shd]=resistance_functions.PrecalculateStomatalResistanceRoof(TempVec_ittm,MeteoData,HumidityAtm,...
            ParVegRoof,SoilPotW_ittm,CiCO2Leaf_ittm,PropOpticalRoof,ra,rb);
    else
        [rs_sun,rs_shd,Ci_sun,Ci_shd]=resistance_functions.PrecalculateStomatalResistanceRoof(TempVec_ittm,MeteoData,HumidityAtm,...
        ParVegRoof,SoilPotW_ittm,CiCO2Leaf_ittm,PropOpticalRoof,RES.raRooftoAtm(ittn-1,1,ittm),RES.rb_LRoof(ittn-1,1,ittm));
    end
else
    rs_sun=Inf;  rs_shd=Inf; Ci_sun=0; Ci_shd=0;
end

rsRoofPreCalc.rs_sun = rs_sun;
rsRoofPreCalc.rs_shd = rs_shd;
rsRoofPreCalc.Ci_sun = Ci_sun;
rsRoofPreCalc.Ci_shd = Ci_shd;

% Precalculate stomatal resistance for tree and for ground vegetation in
% canyon
if ittn==1
    rb_LGround = 50; rb_HGround = 50; rap_can = 100; rap_Htree_In = 100;
    
    [rs_sun_H,rs_shd_H,Ci_sun_H,Ci_shd_H,rs_sun_L,rs_shd_L,Ci_sun_L,Ci_shd_L]=resistance_functions.PrecalculateStomatalResistanceGroundTree(...
             TempVec_ittm,Humidity_ittm,ParVegGround,SoilPotW_ittm,CiCO2Leaf_ittm,...
             MeteoData,geometry,FractionsGround,ParTree,...
		     PropOpticalGround,PropOpticalWall,PropOpticalTree,ParVegTree,...
             SunPosition,ViewFactor,ParWindows,BEM_on,...
             rb_LGround,rb_HGround,rap_can,rap_Htree_In);
else
    [rs_sun_H,rs_shd_H,Ci_sun_H,Ci_shd_H,rs_sun_L,rs_shd_L,Ci_sun_L,Ci_shd_L]=resistance_functions.PrecalculateStomatalResistanceGroundTree(...
         TempVec_ittm,Humidity_ittm,ParVegGround,SoilPotW_ittm,CiCO2Leaf_ittm,...
         MeteoData,geometry,FractionsGround,ParTree,...
	     PropOpticalGround,PropOpticalWall,PropOpticalTree,ParVegTree,...
         SunPosition,ViewFactor,ParWindows,BEM_on,...
         RES.rb_LGround(ittn-1,1,ittm),RES.rb_HGround(ittn-1,1,ittm),RES.rap_can(ittn-1,1,ittm), RES.rap_Htree_In(ittn-1,1,ittm));
end

rsGroundPreCalc.rs_sun_L = rs_sun_L;
rsGroundPreCalc.rs_shd_L = rs_shd_L;
rsGroundPreCalc.Ci_sun_L = Ci_sun_L;
rsGroundPreCalc.Ci_shd_L = Ci_shd_L;

rsTreePreCalc.rs_sun_H = rs_sun_H;
rsTreePreCalc.rs_shd_H = rs_shd_H;
rsTreePreCalc.Ci_sun_H = Ci_sun_H;
rsTreePreCalc.Ci_shd_H = Ci_shd_H;




