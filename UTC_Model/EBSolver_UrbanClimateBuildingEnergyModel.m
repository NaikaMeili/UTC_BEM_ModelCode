function[Ytot]=EBSolver_UrbanClimateBuildingEnergyModel(TemperatureTot,...
        TempVec_ittm,TempVecB_ittm,Humidity_ittm,MeteoData,...
		Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,CiCO2Leaf_ittm,TempDamp_ittm,ViewFactor,...
		Gemeotry_m,ParTree,geometry,FractionsGround,FractionsRoof,...
		WallLayers,ParSoilGround,ParInterceptionTree,...
		PropOpticalGround,PropOpticalWall,PropOpticalTree,...
		ParThermalGround,ParThermalWall,ParVegGround,ParVegTree,...
        ParSoilRoof,PropOpticalRoof,ParThermalRoof,ParVegRoof,...
		SunPosition,HumidityAtm,Anthropogenic,ParCalculation,...
        PropOpticalIndoors,ParHVAC,ParThermalBulidFloor,ParWindows,BEM_on,...
        RESPreCalc,fconvPreCalc,fconv,rsRoofPreCalc,rsGroundPreCalc,rsTreePreCalc,HVACSchedule)


% Temperature vector:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TemperatureTot(:,1)    =   Troof_imp
% TemperatureTot(:,2)    =   Troof_veg
% TemperatureTot(:,3)    = Troof_interior_imp
% TemperatureTot(:,4)    = Troof_interior_veg
%--------------------------------------------------------------
% TemperatureTot(:,5)    =	Temperature ground impervious area
% TemperatureTot(:,6)    =	Temperature ground bare area
% TemperatureTot(:,7)    =	Temperature ground vegetated area
% TemperatureTot(:,8)    =	Temperature sunlit area
% TemperatureTot(:,9)    =	Temperature shaded area
% TemperatureTot(:,10)    =	Temperature tree canopy
% TemperatureTot(:,11)   =	Interior temperature sunlit wall
% TemperatureTot(:,12)   =	Interior temperature shaded wall
% TemperatureTot(:,13)    =	Temperature canyon
% TemperatureTot(:,14)	=	specific humidity canyon
%--------------------------------------------------------------
% TemperatureTot(:,15)   =	Temperature ceiling
% TemperatureTot(:,16)   =	Temperature sunlit wall
% TemperatureTot(:,17)   =	Temperature shaded wall
% TemperatureTot(:,18)   =	Temperature windows
% TemperatureTot(:,19)   =	Temperature ground
% TemperatureTot(:,20)   =	Temperature building internal mass
% TemperatureTot(:,21)   =	Temperature air
% TemperatureTot(:,22)   =	Humidity air
%--------------------------------------------------------------

%--------------------------------------------------------------------------
% TemperatureR(1,1) = Troof_imp
% TemperatureR(1,2) = Troof_veg
% TemperatureR(1,3) = Troof_interior_imp
% TemperatureR(1,4) = Troof_interior_veg

TemperatureR(:,1) = TemperatureTot(:,1);
TemperatureR(:,2) = TemperatureTot(:,2);
TemperatureR(:,3) = TemperatureTot(:,3);
TemperatureR(:,4) = TemperatureTot(:,4);

%--------------------------------------------------------------------------
% TemperatureC(:,1)		=	Temperature ground impervious area
% TemperatureC(:,2)		=	Temperature ground bare area
% TemperatureC(:,3)		=	Temperature ground vegetated area
% TemperatureC(:,4)		=	Temperature sunlit area
% TemperatureC(:,5)		=	Temperature shaded area
% TemperatureC(:,6)		=	Temperature tree canopy
% TemperatureC(:,7)		=	Interior temperature sunlit wall
% TemperatureC(:,8)		=	Interior temperature shaded wall
% TemperatureC(:,9)		=	Temperature canyon
% TemperatureC(:,10)	=	specific humidity canyon

TemperatureC(:,1)		=	TemperatureTot(:,5);
TemperatureC(:,2)		=	TemperatureTot(:,6);
TemperatureC(:,3)		=	TemperatureTot(:,7);
TemperatureC(:,4)		=	TemperatureTot(:,8);
TemperatureC(:,5)		=	TemperatureTot(:,9);
TemperatureC(:,6)		=	TemperatureTot(:,10);
TemperatureC(:,7)		=	TemperatureTot(:,11);
TemperatureC(:,8)		=	TemperatureTot(:,12);
TemperatureC(:,9)		=	TemperatureTot(:,13);
TemperatureC(:,10)	    =	TemperatureTot(:,14);

%--------------------------------------------------------------------------
% TemperatureB(:,1)		=	Temperature ceiling
% TemperatureB(:,2)		=	Temperature sunlit wall
% TemperatureB(:,3)		=	Temperature shaded wall
% TemperatureB(:,4)		=	Temperature window
% TemperatureB(:,5)		=	Temperature ground
% TemperatureB(:,6)		=	Temperature internal mass
% TemperatureB(:,7)		=	Temperature air
% TemperatureB(:,8)		=	Humidity air

TemperatureB(:,1)		=	TemperatureTot(:,15);
TemperatureB(:,2)		=	TemperatureTot(:,16);
TemperatureB(:,3)		=	TemperatureTot(:,17);
TemperatureB(:,4)		=	TemperatureTot(:,18);
TemperatureB(:,5)		=	TemperatureTot(:,19);
TemperatureB(:,6)		=	TemperatureTot(:,20);
TemperatureB(:,7)		=	TemperatureTot(:,21);
TemperatureB(:,8)		=	TemperatureTot(:,22);



[Yroof,G2Roof]=EBSolver_roof(TemperatureR,TemperatureB,TempVec_ittm,MeteoData,...
				Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,CiCO2Leaf_ittm,...
				Gemeotry_m,FractionsRoof,ParSoilRoof,PropOpticalRoof,ParThermalRoof,ParVegRoof,...
				HumidityAtm,Anthropogenic,ParCalculation,BEM_on,RESPreCalc,rsRoofPreCalc);

[Ycanyon,G2WallSun,G2WallShade,SWRabs_t]=EBSolver_canyon(TemperatureC,TemperatureB,TempVec_ittm,Humidity_ittm,MeteoData,...
		Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,CiCO2Leaf_ittm,TempDamp_ittm,ViewFactor,...
		Gemeotry_m,ParTree,geometry,FractionsGround,...
		WallLayers,ParSoilGround,ParInterceptionTree,...
		PropOpticalGround,PropOpticalWall,PropOpticalTree,...
		ParThermalGround,ParThermalWall,ParVegGround,ParVegTree,...
		SunPosition,HumidityAtm,Anthropogenic,ParCalculation,...
        TempVecB_ittm,G2Roof,PropOpticalIndoors,ParHVAC,ParThermalBulidFloor,...
        ParWindows,BEM_on,RESPreCalc,fconvPreCalc,fconv,rsGroundPreCalc,rsTreePreCalc,HVACSchedule);

SWRinWsun = SWRabs_t.SWRabsWallSunTransmitted;
SWRinWshd = SWRabs_t.SWRabsWallShadeTransmitted;

[YBuildInt]=BuildingEnergyModel.EBSolver_Building(TemperatureC,TemperatureB,TempVecB_ittm,TempVec_ittm,Humidity_ittm,MeteoData,...
    SWRinWsun,SWRinWshd,G2Roof,G2WallSun,G2WallShade,TempDamp_ittm,SWRabs_t,...
    Gemeotry_m,PropOpticalIndoors,ParHVAC,ParCalculation,ParThermalBulidFloor,ParWindows,BEM_on,HVACSchedule);

Ytot = [Yroof, Ycanyon, YBuildInt];


