NameOutput	=	'SGTest';
load(['Calculation', NameOutput,'.mat'])


DateTime		=	MeteoDataRaw.Date(1:n);
TempData.Time	=	DateTime(1);
TempData.Zatm	=	Zatm;

for ittm=1:m

[Gemeotry_m,ParTree,geometry,FractionsRoof,FractionsGround,...
	WallLayers,ParSoilRoof,ParSoilGround,ParInterceptionTree,...
	PropOpticalRoof,PropOpticalGround,PropOpticalWall,PropOpticalTree,...
	ParThermalRoof,ParThermalGround,ParThermalWall,ParThermalTree,...
	ParVegRoof,ParVegGround,ParVegTree]=feval(strcat('data_functions.Data_UEHM_site_',Name_Site),TempData,ittm,NaN);

clearvars -except geometry FractionsGround FractionsRoof m ittm MeteoDataRaw MeteoDataPreston_h DateTime NameOutput
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['Calculation', NameOutput,'.mat'])

%% Energy blance by surface including solver to check for assignment errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EB_RoofImpSolv		=	SWRabs.SWRabsRoofImp(:,:,ittm) + LWRabs.LWRabsRoofImp(:,:,ittm) - Gflux.G1RoofImp(:,:,ittm) - Hflux.HfluxRoofImp(:,:,ittm) - LEflux.LEfluxRoofImp(:,:,ittm) - Solver.ValuesEB(:,1,ittm);
EB_RoofVegSolv		=	SWRabs.SWRabsRoofVeg(:,:,ittm) + LWRabs.LWRabsRoofVeg(:,:,ittm) - Gflux.G1RoofVeg(:,:,ittm) - Hflux.HfluxRoofVeg(:,:,ittm) - LEflux.LEfluxRoofVeg(:,:,ittm) - Solver.ValuesEB(:,2,ittm);

EB_GroundImpSolv	=	SWRabs.SWRabsGroundImp(:,:,ittm) + LWRabs.LWRabsGroundImp(:,:,ittm) - Gflux.G1GroundImp(:,:,ittm) - Hflux.HfluxGroundImp(:,:,ittm) - LEflux.LEfluxGroundImp(:,:,ittm) - Solver.ValuesEB(:,3,ittm);
EB_GroundBareSolv	=	SWRabs.SWRabsGroundBare(:,:,ittm) + LWRabs.LWRabsGroundBare(:,:,ittm) - Gflux.G1GroundBare(:,:,ittm) - Hflux.HfluxGroundBare(:,:,ittm) - LEflux.LEfluxGroundBare(:,:,ittm) - Solver.ValuesEB(:,4,ittm);
EB_GroundVegSolv	=	SWRabs.SWRabsGroundVeg(:,:,ittm) + LWRabs.LWRabsGroundVeg(:,:,ittm) - Gflux.G1GroundVeg(:,:,ittm) - Hflux.HfluxGroundVeg(:,:,ittm) - LEflux.LEfluxGroundVeg(:,:,ittm) - Solver.ValuesEB(:,5,ittm);

EB_TreeSolv			=	SWRabs.SWRabsTree(:,:,ittm) + LWRabs.LWRabsTree(:,:,ittm) - Hflux.HfluxTree(:,:,ittm) - LEflux.LEfluxTree(:,:,ittm) - Solver.ValuesEB(:,8,ittm);
EB_WallSunSolv		=	SWRabs.SWRabsWallSun(:,:,ittm) + LWRabs.LWRabsWallSun(:,:,ittm) - Gflux.G1WallSun(:,:,ittm) - Hflux.HfluxWallSun(:,:,ittm) - LEflux.LEfluxWallSun(:,:,ittm) - Solver.ValuesEB(:,6,ittm);
EB_WallShadeSolv	=	SWRabs.SWRabsWallShade(:,:,ittm) + LWRabs.LWRabsWallShade(:,:,ittm) - Gflux.G1WallShade(:,:,ittm) - Hflux.HfluxWallShade(:,:,ittm) - LEflux.LEfluxWallShade(:,:,ittm) - Solver.ValuesEB(:,7,ittm);

EB_HfluxCanyonSolv	=	FractionsGround.fimp.*Hflux.HfluxGroundImp(:,:,ittm) + FractionsGround.fbare.*Hflux.HfluxGroundBare(:,:,ittm) + FractionsGround.fveg.*Hflux.HfluxGroundVeg(:,:,ittm) ...
						+ 4.*geometry.radius_tree.*Hflux.HfluxTree(:,:,ittm) + geometry.hcanyon.*Hflux.HfluxWallSun(:,:,ittm) + geometry.hcanyon.*Hflux.HfluxWallShade(:,:,ittm) + Anthropo.Qf_canyon(:,:,ittm)...
						- Hflux.HfluxCanyon(:,:,ittm) - Solver.ValuesEB(:,9,ittm);
					
EB_LEfluxCanyonSolv	=	FractionsGround.fimp.*LEflux.LEfluxGroundImp(:,:,ittm) + FractionsGround.fbare.*LEflux.LEfluxGroundBare(:,:,ittm) + FractionsGround.fveg.*LEflux.LEfluxGroundVeg(:,:,ittm) ...
						+ 4.*geometry.radius_tree.*LEflux.LEfluxTree(:,:,ittm) + geometry.hcanyon.*LEflux.LEfluxWallSun(:,:,ittm) + geometry.hcanyon.*LEflux.LEfluxWallShade(:,:,ittm) - LEflux.LEfluxCanyon(:,:,ittm) - Solver.ValuesEB(:,10,ittm);
				
MaxEB_RoofImpSolv		=	max(abs(EB_RoofImpSolv));
MaxEB_RoofVegSolv		=	max(abs(EB_RoofVegSolv));
MaxEB_GroundImpSolv		=	max(max(abs(EB_GroundImpSolv)));
MaxEB_GroundBareSolv	=	max(max(abs(EB_GroundBareSolv)));
MaxEB_GroundVegSolv		=	max(max(abs(EB_GroundVegSolv)));
MaxEB_TreeSolv			=	max(max(abs(EB_TreeSolv)));
MaxEB_WallSunSolv		=	max(max(abs(EB_WallSunSolv)));
MaxEB_WallShadeSolv		=	max(max(abs(EB_WallShadeSolv)));
MaxEB_HfluxCanyonSolv	=	max(max(abs(EB_HfluxCanyonSolv)));
MaxEB_LEfluxCanyonSolv	=	max(max(abs(EB_LEfluxCanyonSolv)));

MaxEB_Solv	=	max([MaxEB_RoofImpSolv,MaxEB_RoofVegSolv,MaxEB_GroundImpSolv,...
				MaxEB_GroundBareSolv,MaxEB_GroundVegSolv,MaxEB_TreeSolv,...
				MaxEB_WallSunSolv,MaxEB_WallShadeSolv,MaxEB_HfluxCanyonSolv,MaxEB_LEfluxCanyonSolv]);

if any(MaxEB_Solv>1e-6)
	disp(['Please check surface energy balance. There is an assignment error. For further details see EnergyBalanceCheck.m at ittm ' num2str(ittm)])
end

%% Check solver convergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SolverRoofSum			=	sum(Solver.ValuesEB(:,1:2,ittm),1);
SolverRoofSumAbs		=	sum(abs(Solver.ValuesEB(:,1:2,ittm)),1);
SolverRoofMaxAbs		=	max(abs(Solver.ValuesEB(:,1:2,ittm)),[],1);
SolverRoofConvCount		=	sum(abs(Solver.ValuesEB(:,1:2,ittm))>10^-6,1);

SolverCanyonSum			=	sum(Solver.ValuesEB(:,3:10,ittm),1);
SolverCanyonSumAbs		=	sum(abs(Solver.ValuesEB(:,3:10,ittm)),1);
SolverCanyonMaxAbs		=	max(abs(Solver.ValuesEB(:,3:10,ittm)),[],1);
SolverCanyonConvCount	=	sum(abs(Solver.ValuesEB(:,3:10,ittm))>10^-6,1);

if any(SolverCanyonConvCount>0)
	disp(['The solver has problems converging at ittm = ',num2str(ittm)])
	disp(['The overall sum of the energy balance difference of each canyon component is ',num2str(SolverCanyonSum),' W/m^2.'])
	disp(['The mean energy balance difference of each canyon component is ',num2str(SolverCanyonSum./length(Solver.ValuesEB(:,3:10,ittm))),' W/m^2.'])
	disp(['The number of nonconvergence of each canyon component is ',num2str(SolverCanyonConvCount)])
	
	itteration	=	length(Solver.Success(:,1,m));
	colors		=	[cellstr('r');cellstr('r:');cellstr('b');cellstr('b:');cellstr('b.');...
				cellstr('g');cellstr('m');cellstr('c');cellstr('m:');cellstr('c:');...
				cellstr('k');cellstr('y')];
			
	EBNames	=	fieldnames(EB);

	figure
	for i=1:size(EBNames,1)
		data	=	EB.(cell2mat(EBNames(i)))(:,1,ittm);
		name	=	cell2mat(EBNames(i));
		plot(DateTime,data,colors{i},'DisplayName',name)
		hold on	
	end
	legend('show')
	xlim([DateTime(1,1) DateTime(itteration,1)])

end

if any(SolverRoofConvCount>0)
	disp(['The solver has problems converging at ittm = ',num2str(ittm)])
	disp(['The overall sum of the energy balance difference of each roof component is ',num2str(SolverRoofSum),' W/m^2.'])
	disp(['The mean energy balance difference of each roof component is ',num2str(SolverRoofSum./length(Solver.ValuesEB(:,1:2,ittm))),' W/m^2.'])
	disp(['The number of nonconvergence of each roof component is ',num2str(SolverRoofConvCount)])
	
	itteration	=	length(Solver.Success(:,1,m));
	%DateTime	=	datetime(MeteoDataRaw.T_atm(1:itteration,1),1,MeteoDataRaw.T_atm(1:itteration,2),MeteoDataRaw.T_atm(1:itteration,3), 0, 0);
	colors		=	[cellstr('r');cellstr('r:');cellstr('b');cellstr('b:');cellstr('b.');...
				cellstr('g');cellstr('m');cellstr('c');cellstr('m:');cellstr('c:');...
				cellstr('k');cellstr('y')];
			
	EBNames	=	fieldnames(EB);

	figure
	for i=1:size(EBNames,1)
		data	=	EB.(cell2mat(EBNames(i)))(:,1,ittm);
		name	=	cell2mat(EBNames(i));
		plot(DateTime,data,colors{i},'DisplayName',name)
		hold on	
	end
	legend('show')
	xlim([DateTime(1,1) DateTime(itteration,1)])
	
end

SolverConvergenceCan(:,:,ittm)	=	FractionsGround.fimp.*Solver.ValuesCanyon(:,1,ittm) + FractionsGround.fbare.*Solver.ValuesCanyon(:,2,ittm) + FractionsGround.fveg.*Solver.ValuesCanyon(:,3,ittm) ...
							+ 4.*geometry.radius_tree.*Solver.ValuesCanyon(:,6,ittm) + geometry.hcanyon.*Solver.ValuesCanyon(:,4,ittm) + geometry.hcanyon.*Solver.ValuesCanyon(:,5,ittm) ...
							+ Solver.ValuesCanyon(:,9,ittm) + Solver.ValuesCanyon(:,10,ittm);
						
SolverCanyonMax			=	max(SolverConvergenceCan(:,:,ittm),[],1);
SolverCanyonMean		=	nanmean(SolverConvergenceCan(:,:,ittm),1);

if any(SolverCanyonMax>10^-3)
	disp(['The maximum canyon solver convergence issue is ',num2str(SolverCanyonMax),' W/m^2.'])
	disp(['The mean canyon solver convergence is ',num2str(SolverCanyonMean),' W/m^2.'])
end
	
% Overall energy balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ShortwaveRadiationAbsRoof=	FractionsRoof.fimp.*SWRabs.SWRabsRoofImp(:,:,ittm) + FractionsRoof.fveg.*SWRabs.SWRabsRoofVeg(:,:,ittm);
LongwaveRadiationRoof	=	FractionsRoof.fimp.*LWRabs.LWRabsRoofImp(:,:,ittm) + FractionsRoof.fveg.*LWRabs.LWRabsRoofVeg(:,:,ittm);						
ConductiveHeatG1Roof	=	FractionsRoof.fimp.*Gflux.G1RoofImp(:,:,ittm) + FractionsRoof.fveg.*Gflux.G1RoofVeg(:,:,ittm);					
SensibleHeatRoof		=	FractionsRoof.fimp.*Hflux.HfluxRoofImp(:,:,ittm) + FractionsRoof.fveg.*Hflux.HfluxRoofVeg(:,:,ittm);				
LatentHeatRoof			=	FractionsRoof.fimp.*LEflux.LEfluxRoofImp(:,:,ittm) + FractionsRoof.fveg.*LEflux.LEfluxRoofVeg(:,:,ittm);
												
ShortwaveRadiationAbsCan=	FractionsGround.fimp.*SWRabs.SWRabsGroundImp(:,:,ittm) + FractionsGround.fbare.*SWRabs.SWRabsGroundBare(:,:,ittm) + FractionsGround.fveg.*SWRabs.SWRabsGroundVeg(:,:,ittm)...
							+ 4.*geometry.radius_tree.*SWRabs.SWRabsTree(:,:,ittm) + geometry.hcanyon.*SWRabs.SWRabsWallSun(:,:,ittm) + geometry.hcanyon.* SWRabs.SWRabsWallShade(:,:,ittm);

LongwaveRadiationCan	=	FractionsGround.fimp.*LWRabs.LWRabsGroundImp(:,:,ittm) + FractionsGround.fbare.*LWRabs.LWRabsGroundBare(:,:,ittm) + FractionsGround.fveg.*LWRabs.LWRabsGroundVeg(:,:,ittm)...
							+ 4.*geometry.radius_tree.*LWRabs.LWRabsTree(:,:,ittm) + geometry.hcanyon.*LWRabs.LWRabsWallSun(:,:,ittm) + geometry.hcanyon.* LWRabs.LWRabsWallShade(:,:,ittm);
						
ConductiveHeatG1Can		=	FractionsGround.fimp.*Gflux.G1GroundImp(:,:,ittm) + FractionsGround.fbare.*Gflux.G1GroundBare(:,:,ittm) + FractionsGround.fveg.*Gflux.G1GroundVeg(:,:,ittm)...
							+ geometry.hcanyon.*Gflux.G1WallSun(:,:,ittm) + geometry.hcanyon.* Gflux.G1WallShade(:,:,ittm);
						
SensibleHeatCan			=	FractionsGround.fimp.*Hflux.HfluxGroundImp(:,:,ittm) + FractionsGround.fbare.*Hflux.HfluxGroundBare(:,:,ittm) + FractionsGround.fveg.*Hflux.HfluxGroundVeg(:,:,ittm)...
							+ 4.*geometry.radius_tree.*Hflux.HfluxTree(:,:,ittm) + geometry.hcanyon.*Hflux.HfluxWallSun(:,:,ittm) + geometry.hcanyon.*Hflux.HfluxWallShade(:,:,ittm)+ Anthropo.Qf_canyon(:,:,ittm);
				
LatentHeatCan			=	FractionsGround.fimp.*LEflux.LEfluxGroundImp(:,:,ittm) + FractionsGround.fbare.*LEflux.LEfluxGroundBare(:,:,ittm) + FractionsGround.fveg.*LEflux.LEfluxGroundVeg(:,:,ittm) ...
							+ 4.*geometry.radius_tree.*LEflux.LEfluxTree(:,:,ittm) + geometry.hcanyon.*LEflux.LEfluxWallSun(:,:,ittm) + geometry.hcanyon.*LEflux.LEfluxWallShade(:,:,ittm);
					
SWRRoofCalcDif	=	max(ShortwaveRadiationAbsRoof-SWRabs.SWRabsTotalRoof(:,:,ittm));
LWRRoofCalcDif	=	max(LongwaveRadiationRoof-LWRabs.LWRabsTotalRoof(:,:,ittm));
G1RoofCalcDif	=	max(ConductiveHeatG1Roof-Gflux.G1Roof(:,:,ittm));
HRoofCalcDif	=	max(SensibleHeatRoof-Hflux.HfluxRoof(:,:,ittm));
LERoofCalcDif	=	max(LatentHeatRoof-LEflux.LEfluxRoof(:,:,ittm));
SWRCanyonCalcDif=	max(ShortwaveRadiationAbsCan-SWRabs.SWRabsTotalCanyon(:,:,ittm));
LWRCanyonCalcDif=	max(LongwaveRadiationCan-LWRabs.LWRabsTotalCanyon(:,:,ittm));
G1CanyonCalcDif	=	max(ConductiveHeatG1Can-Gflux.G1Canyon(:,:,ittm));
HCanyonCalcDif	=	max(SensibleHeatCan-Hflux.HfluxCanyon(:,:,ittm)- Solver.ValuesCanyon(:,9,ittm));
LECanyonCalcDif	=	max(LatentHeatCan-LEflux.LEfluxCanyon(:,:,ittm)- Solver.ValuesCanyon(:,10,ittm));

CalcDiff	=	max([SWRRoofCalcDif,LWRRoofCalcDif,G1RoofCalcDif,HRoofCalcDif,...
				LERoofCalcDif,SWRCanyonCalcDif,LWRCanyonCalcDif,G1CanyonCalcDif,...
				HCanyonCalcDif,LECanyonCalcDif]);

if any(CalcDiff>1e-6)
	  disp(['There is a variable assignment error. Please check the overall energy balance calculations.'])
end
												
EBRoof1					=	ShortwaveRadiationAbsRoof + LongwaveRadiationRoof - ConductiveHeatG1Roof - SensibleHeatRoof - LatentHeatRoof;
EBRoof2					=	SWRabs.SWRabsTotalRoof(:,:,ittm) + LWRabs.LWRabsTotalRoof(:,:,ittm) - Gflux.G1Roof(:,:,ittm) - Hflux.HfluxRoof(:,:,ittm) - LEflux.LEfluxRoof(:,:,ittm);
EBCanyon1				=	ShortwaveRadiationAbsCan + LongwaveRadiationCan + Anthropo.Qf_canyon(:,:,ittm) - ConductiveHeatG1Can - (SensibleHeatCan - Solver.ValuesCanyon(:,9,ittm)) - (LatentHeatCan - Solver.ValuesCanyon(:,10,ittm)) - SolverConvergenceCan(:,:,ittm);
EBCanyon2				=	SWRabs.SWRabsTotalCanyon(:,:,ittm) + LWRabs.LWRabsTotalCanyon(:,:,ittm) + Anthropo.Qf_canyon(:,:,ittm) - Gflux.G1Canyon(:,:,ittm) - Hflux.HfluxCanyon(:,:,ittm) - LEflux.LEfluxCanyon(:,:,ittm) - SolverConvergenceCan(:,:,ittm);
EBUrban					=	SWRabs.SWRabsTotalUrban(:,:,ittm) + LWRabs.LWRabsTotalUrban(:,:,ittm) + geometry.wcanyon_norm.*Anthropo.Qf_canyon(:,:,ittm) - Gflux.G1Urban(:,:,ittm) - Hflux.HfluxUrban(:,:,ittm) - LEflux.LEfluxUrban(:,:,ittm)- geometry.wcanyon_norm.*SolverConvergenceCan(:,:,ittm);

MaxEBRoof1				=	max(abs(EBRoof1));
MaxEBRoof2				=	max(abs(EBRoof2)); % There is a slight difference. Check why
MaxEBCanyon1			=	max(abs(EBCanyon1));
MaxEBCanyon2			=	max(abs(EBCanyon2)); % There is a slight difference. Check why
MaxEBUrban				=	max(abs(EBUrban));

MaxEBRoof	=	max([MaxEBRoof1,MaxEBRoof2]);
MaxEBCanyon	=	max([MaxEBCanyon1,MaxEBCanyon2]);
MaxEBUrban	=	max([MaxEBUrban]);

if any(MaxEBRoof>10^-3)
	disp(['The maximum energy balance of the roof is ', num2str(MaxEBRoof),' W/m^2 at ittm ', num2str(ittm)])
end

if any(MaxEBCanyon>10^-3)
	disp(['The maximum energy balance of the canyon is ', num2str(MaxEBCanyon),' W/m^2 at ittm ', num2str(ittm)])
end

if any(MaxEBUrban>10^-3)
	disp(['The maximum energy balance of the urban is ', num2str(MaxEBUrban),' W/m^2 at ittm ', num2str(ittm)])
end

% Mean
MeanEBRoof1				=	nanmean(ShortwaveRadiationAbsRoof) + nanmean(LongwaveRadiationRoof) - nanmean(ConductiveHeatG1Roof) - nanmean(SensibleHeatRoof) - nanmean(LatentHeatRoof);
MeanEBRoof2				=	nanmean(SWRabs.SWRabsTotalRoof(:,:,ittm)) + nanmean(LWRabs.LWRabsTotalRoof(:,:,ittm)) - nanmean(Gflux.G1Roof(:,:,ittm)) - nanmean(Hflux.HfluxRoof(:,:,ittm)) - nanmean(LEflux.LEfluxRoof(:,:,ittm));
MeanEBCanyon1			=	nanmean(ShortwaveRadiationAbsCan) + nanmean(LongwaveRadiationCan) + nanmean(Anthropo.Qf_canyon(:,:,ittm)) - nanmean(ConductiveHeatG1Can) - nanmean(SensibleHeatCan - Solver.ValuesCanyon(:,9,ittm)) - nanmean(LatentHeatCan - Solver.ValuesCanyon(:,10,ittm)) - nanmean(SolverConvergenceCan(:,:,ittm));
MeanEBCanyon2			=	nanmean(SWRabs.SWRabsTotalCanyon(:,:,ittm)) + nanmean(LWRabs.LWRabsTotalCanyon(:,:,ittm)) + nanmean(Anthropo.Qf_canyon(:,:,ittm)) - nanmean(Gflux.G1Canyon(:,:,ittm)) - nanmean(Hflux.HfluxCanyon(:,:,ittm)) - nanmean(LEflux.LEfluxCanyon(:,:,ittm)) - nanmean(SolverConvergenceCan(:,:,ittm));
MeanEBUrban				=	nanmean(SWRabs.SWRabsTotalUrban(:,:,ittm)) + nanmean(LWRabs.LWRabsTotalUrban(:,:,ittm)) + nanmean(geometry.wcanyon_norm.*Anthropo.Qf_canyon(:,:,ittm)) - nanmean(Gflux.G1Urban(:,:,ittm)) - nanmean(Hflux.HfluxUrban(:,:,ittm)) - nanmean(LEflux.LEfluxUrban(:,:,ittm))- nanmean(geometry.wcanyon_norm.*SolverConvergenceCan(:,:,ittm));

MeanEBRoof	=	max([MeanEBRoof1,MeanEBRoof2]);
MeanEBCanyon=	max([MeanEBCanyon1,MeanEBCanyon2]);
MeanEBUrban	=	max([MeanEBUrban]);

if any(MeanEBRoof>0.1)
	disp(['The mean energy balance of the roof is ', num2str(MeanEBRoof),' W/m^2 at ittm ', num2str(ittm)])
end

if any(MeanEBCanyon>0.1)
	disp(['The mean energy balance of the canyon is ', num2str(MeanEBCanyon),' W/m^2 at ittm ', num2str(ittm)])
end

if any(MeanEBUrban>0.1)
	disp(['The mean energy balance of the urban is ', num2str(MeanEBUrban),' W/m^2 at ittm ', num2str(ittm)])
end


end
