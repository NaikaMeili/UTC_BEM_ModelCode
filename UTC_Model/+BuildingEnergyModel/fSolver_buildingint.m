function[T,fval,exitflag]=fSolver_buildingint(TempVecB,MeteoData,HumidityAtm,SWRinWsun,SWRinWshd,...
    Tairout,qairout,G2Roof,G2WallSun,G2WallShade,TBdamp_ittm,...
    Gemeotry_m,PropOpticalIndoors,ParHVAC,ParCalculation,ParThermalBulidFloor)


% Nonlinear system solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Temperature vector:
% TemperatureB(:,1)		=	Temperature ceiling
% TemperatureB(:,2)		=	Temperature sunlit wall
% TemperatureB(:,3)		=	Temperature shaded wall
% TemperatureB(:,4)		=	Temperature ground
% TemperatureB(:,5)		=	Temperature air
% TemperatureB(:,6)		=	Humidity air


Opt_Solv = optimoptions('lsqnonlin','Display','off','FunctionTolerance',10^-10,'MaxFunctionEvaluations',300);


% Use temperature from previous time step as a starting point
TemperatureB	=	[TempVecB.Tceiling, TempVecB.Tinwallsun,TempVecB.Tinwallshd,...
					TempVecB.Tinground,TempVecB.Tbin,TempVecB.qbin];

if ParHVAC.ACon == 1
    TemperatureB	=	[TempVecB.Tceiling, TempVecB.Tinwallsun,TempVecB.Tinwallshd,...
					TempVecB.Tinground,ParHVAC.Tsetpoint,TempVecB.qbin];
end


DeltaT = 30;
lb	=	[TempVecB.Tceiling-DeltaT, TempVecB.Tinwallsun-DeltaT,TempVecB.Tinwallshd-DeltaT,...
					TempVecB.Tinground-DeltaT,TempVecB.Tbin-DeltaT,0];
ub	=	[TempVecB.Tceiling+DeltaT, TempVecB.Tinwallsun+DeltaT,TempVecB.Tinwallshd+DeltaT,...
					TempVecB.Tinground+DeltaT,TempVecB.Tbin+DeltaT,1];

[T,~,fval,exitflag] = lsqnonlin(@targetFun,TemperatureB,lb,ub,Opt_Solv);

%--------------------------------------------------------------------------
% If solver failed, retry with different starting value
if sum(abs(fval)>0.01)>0
	TT = MeteoData.Tatm;
	TemperatureB	=	[TT,TT,TT,TT,TT,HumidityAtm.AtmSpecific];
    
    [T,~,fval,exitflag] = lsqnonlin(@targetFun,TemperatureB,lb,ub,Opt_Solv);
end

for i=1:3
	if sum(abs(fval)>0.01)>0
        TemperatureB	=	[TempVecB.Tceiling+i, TempVecB.Tinwallsun+i,TempVecB.Tinwallshd+i,...
					TempVecB.Tinground+i,TempVecB.Tbin+i,HumidityAtm.AtmSpecific];

		[T,~,fval,exitflag] = lsqnonlin(@targetFun,TemperatureB,lb,ub,Opt_Solv);
	else
		break
	end	
end

for i=1:3
	if sum(abs(fval)>0.01)>0
		TemperatureB	=	[TempVecB.Tceiling-i, TempVecB.Tinwallsun-i,TempVecB.Tinwallshd-i,...
					TempVecB.Tinground-i,TempVecB.Tbin-i,HumidityAtm.AtmSpecific];

		[T,~,fval,exitflag] = lsqnonlin(@targetFun,TemperatureB,lb,ub,Opt_Solv);
	else
		break
	end	
end

%--------------------------------------------------------------------------
% Create dummy function to pass further variables to the solver
function YBuildInt = targetFun(TemperatureB)
	YBuildInt = BuildingEnergyModel.EBSolver_BuildingInt(TemperatureB,MeteoData,SWRinWsun,SWRinWshd,...
    TempVecB,Tairout,qairout,G2Roof,G2WallSun,G2WallShade,TBdamp_ittm,...
    Gemeotry_m,PropOpticalIndoors,ParHVAC,ParCalculation,ParThermalBulidFloor);
end

end


