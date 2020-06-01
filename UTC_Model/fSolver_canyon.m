function[T,fval,exitflag]=fSolver_canyon(TempVec,Humidity,MeteoData,itt,...
		Int,ExWater,Vwater,Owater,SoilPotW,CiCO2Leaf,TempDamp,ViewFactor,Opt_Solv,...
		Gemeotry_m,ParTree,geometry,FractionsGround,...
		WallLayers,ParSoilGround,ParInterceptionTree,...
		PropOpticalGround,PropOpticalWall,PropOpticalTree,...
		ParThermalGround,ParThermalWall,ParVegGround,ParVegTree,...
		SunPosition,HumidityAtm,Anthropogenic,ParCalculation)

% Nonlinear system solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temperature vector:
% TempVec(:,1)	=	Temperature ground impervious area
% TempVec(:,2)	=	Temperature ground bare area
% TempVec(:,3)	=	Temperature ground vegetated area
% TempVec(:,4)	=	Temperature sunlit wall
% TempVec(:,5)	=	Temperature shaded wall
% TempVec(:,6)	=	Temperature tree canopy
% TempVec(:,7)	=	Interior temperature sunlit wall
% TempVec(:,8)	=	Interior temperature shaded wall
% TempVec(:,9)	=	Temperature canyon air
% TempVec(:,10)	=	Humidity canyon air


Opt_Solv = optimoptions('lsqnonlin','Display','off','FunctionTolerance',10^-10,'MaxFunctionEvaluations',300);


% Use temperature from previous time step as a starting point
TemperatureC	=	[TempVec.TGroundImp(itt,1), TempVec.TGroundBare(itt,1),TempVec.TGroundVeg(itt,1),...
					TempVec.TWallSun(itt,1),TempVec.TWallShade(itt,1),TempVec.TTree(itt,1),TempVec.TWallIntSun(itt,1),...
					TempVec.TWallIntShade(itt,1),TempVec.TCanyon(itt,1),Humidity.CanyonSpecific(itt,1)];
							
lb	=	[243,243,243,243,243,243,243,243,243,0];
% ub	=	[inf,inf,inf,inf,inf,inf,inf,inf,inf,inf];

[T,~,fval,exitflag] = lsqnonlin(@targetFun,TemperatureC,lb,[],Opt_Solv);

% If solver failed, retry with different starting value
if sum(abs(fval)>0.01)>0
	TT = MeteoData.Tatm(itt);
	TemperatureC	=	[TT,TT,TT,TT,TT,TT,TT,TT,TT,Humidity.AtmSpecific(itt)];
    
    [T,~,fval,exitflag] = lsqnonlin(@targetFun,TemperatureC,lb,[],Opt_Solv);
end

for i=1:3
	if sum(abs(fval)>0.01)>0
		TemperatureC	=	[TempVec.TGroundImp(itt,1)+i, TempVec.TGroundBare(itt,1)+i,TempVec.TGroundVeg(itt,1)+i,...
						TempVec.TWallSun(itt,1)+i,TempVec.TWallShade(itt,1)+i,TempVec.TTree(itt,1)+i,TempVec.TWallIntSun(itt,1)+i,...
						TempVec.TWallIntShade(itt,1)+i,TempVec.TCanyon(itt,1)+i,Humidity.CanyonSpecific(itt,1)];

		[T,~,fval,exitflag] = lsqnonlin(@targetFun,TemperatureC,lb,[],Opt_Solv);
	else
		break
	end	
end

for i=1:3
	if sum(abs(fval)>0.01)>0
		TemperatureC	=	[TempVec.TGroundImp(itt,1)-i, TempVec.TGroundBare(itt,1)-i,TempVec.TGroundVeg(itt,1)-i,...
						TempVec.TWallSun(itt,1)-i,TempVec.TWallShade(itt,1)-i,TempVec.TTree(itt,1)-i,TempVec.TWallIntSun(itt,1)-i,...
						TempVec.TWallIntShade(itt,1)-i,TempVec.TCanyon(itt,1)-i,Humidity.CanyonSpecific(itt,1)];

		[T,~,fval,exitflag] = lsqnonlin(@targetFun,TemperatureC,lb,[],Opt_Solv);
	else
		break
	end	
end


% Create dummy function to pass further variables to the solver
function Ycanyon = targetFun(TemperatureC)
	Ycanyon = EBSolver_canyon(TemperatureC,TempVec,MeteoData,itt,...
	Int,ExWater,Vwater,Owater,SoilPotW,CiCO2Leaf,TempDamp,ViewFactor,...
	Gemeotry_m,ParTree,geometry,FractionsGround,...
	WallLayers,ParSoilGround,ParInterceptionTree,...
	PropOpticalGround,PropOpticalWall,PropOpticalTree,...
	ParThermalGround,ParThermalWall,ParVegGround,ParVegTree,...
	SunPosition,HumidityAtm,Anthropogenic,ParCalculation);
end

end
