function[T,fval,exitflag]=fSolver_roof(TempVec,MeteoData,itt,...
		Int,ExWater,Vwater,Owater,SoilPotW,CiCO2Leaf,Opt_Solv,...
		Gemeotry_m,FractionsRoof,ParSoilRoof,PropOpticalRoof,ParThermalRoof,ParVegRoof,...
		HumidityAtm,Anthropogenic,ParCalculation)

% Nonlinear system solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temperature vector:
% TempVec(:,1)	=	Temperature surface impervious area
% TempVec(:,2)	=	Temperature surface vegetated area
% TempVec(:,3)	=	interior temperature impervious area
% TempVec(:,4)	=	interior temperature vegetated area
% Use temperature from previous time step as a starting point

% Use temperature from previous time step as a starting point
TemperatureR	=	[TempVec.TRoofImp(itt,1), TempVec.TRoofVeg(itt,1),...
					TempVec.TRoofIntImp(itt,1),TempVec.TRoofIntVeg(itt,1)];

lb	=	[243,243,243,243];
% ub	=	[373,373,373,373];

[T,~,fval,exitflag] = lsqnonlin(@targetFun,TemperatureR,lb,[],Opt_Solv);

if exitflag~=1
    TT = MeteoData.Tatm(itt);
    TemperatureC = [TT, TT, TT, TT];
    [T,~,fval,exitflag] = lsqnonlin(@targetFun,TemperatureC,lb,[],Opt_Solv);
end

% Create dummy function to pass further variables to the solver
function Yroof = targetFun(TemperatureR)
		Yroof = EBSolver_roof(TemperatureR,TempVec,MeteoData,itt,...
				Int,ExWater,Vwater,Owater,SoilPotW,CiCO2Leaf,...
				Gemeotry_m,FractionsRoof,ParSoilRoof,PropOpticalRoof,ParThermalRoof,ParVegRoof,...
				HumidityAtm,Anthropogenic,ParCalculation);
	end

end
