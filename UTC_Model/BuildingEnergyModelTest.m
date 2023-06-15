SimC = load('CalculationSG.mat');

Name_SiteFD	=	'SG';

n = length(SimC.MeteoDataRaw.Date);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% Initialize temperature calculation vector
TempVecBNames	=	{'Tceiling';'Tinwallsun';'Tinwallshd';'Tinground';'Tbin';'qbin'};

for i=1:size(TempVecBNames,1)
	TempVecB.(cell2mat(TempVecBNames(i)))			=	zeros(n,1);
	TempVecB.(cell2mat(TempVecBNames(i)))(1,:,:)	=	SimC.MeteoDataRaw.T_atm(1,1,1);
end

TempVecB.qbin(1,1,1) = SimC.Humidity.AtmSpecific(1,1,1);


% Initialize energy flux outputs
HbuildIntNames	=	{'HBinRoof';'HbinWallSun';'HbinWallshd';'HBinGround';...
    'HbuildInSurf';'Hvent';'Hequip';'Hpeople';'H_AC';'H_air'};

for i=1:size(HbuildIntNames,1)
	HbuildInt.(cell2mat(HbuildIntNames(i)))			=	zeros(n,1);
end

LEbuildIntNames	=	{'LEvent';'LEequip';'LEpeople';'LE_AC';'LE_air'};

for i=1:size(LEbuildIntNames,1)
	LEbuildInt.(cell2mat(LEbuildIntNames(i)))			=	zeros(n,1);
end

GbuildIntNames	=	{'G2Roof';'G2WallSun';'G2WallShade';'Gfloor'};

for i=1:size(GbuildIntNames,1)
	GbuildInt.(cell2mat(GbuildIntNames(i)))			=	zeros(n,1);
end

SWRabsBNames	=	{'SWRabsCeiling';'SWRabsWallsun';'SWRabsWallshade';'SWRabsGround'};

for i=1:size(SWRabsBNames,1)
	SWRabsB.(cell2mat(SWRabsBNames(i)))			=	zeros(n,1);
end

LWRabsBNames	=	{'LWRabsCeiling';'LWRabsWallsun';'LWRabsWallshade';'LWRabsGround'};

for i=1:size(LWRabsBNames,1)
	LWRabsB.(cell2mat(LWRabsBNames(i)))			=	zeros(n,1);
end

TBdamp.TingroundDamp		=	zeros(n,1);
TBdamp.TingroundDamp(1,1,1)	=	mean(SimC.MeteoDataRaw.T_atm);
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Things we can keep constant for now
%----------------------------------------------
% Urban geometry (in metres)
Gemeotry_m.Height_canyon  =   10;
Gemeotry_m.Width_roof   =   3;

% albedo/emissivities
PropOpticalIndoors.abc = 0.1;
PropOpticalIndoors.abw = 0.2;
PropOpticalIndoors.abg = 0.3;
PropOpticalIndoors.ec = 0.95;
PropOpticalIndoors.eg = 0.95;
PropOpticalIndoors.ew = 0.95;

ParThermalBulidFloor.lan_dry_imp = 1.2;		% Thermal conductivity dry solid [W/m K]
ParThermalBulidFloor.cv_s_imp = 1.5*10^6;	% Volumetric heat capacity solid [J/m^3 K]

% Aircon parameters
ParHVAC.ACon = 1;
ParHVAC.ACH = 1; % air changes per hour (1/h)
ParHVAC.Tsetpoint = 273.15+25;
ParHVAC.COP = 2.5;

% Temporal resolution
ParCalculation.dth = 1; % (h)
ParCalculation.dts = 3600; % (s) number of seconds in 1 hour


ittm = 1;
for ittn	= 1:n  

    [SunPosition,MeteoData,HumidityAtm,Anthropogenic,location,ParCalculation]...
	=feval(strcat('data_functions.UEHMForcingData_',Name_SiteFD),SimC.MeteoDataRaw,ittn);

    for i=1:size(TempVecBNames,1)
		if ittn==1
			TempVecB_ittm.(cell2mat(TempVecBNames(i)))	=	TempVecB.(cell2mat(TempVecBNames(i)))(1,:,ittm); 
		else
			TempVecB_ittm.(cell2mat(TempVecBNames(i)))	=	TempVecB.(cell2mat(TempVecBNames(i)))(ittn-1,:,ittm);
		end
    end

    if ittn==1
        TBdamp_ittm.TingroundDamp = TBdamp.TingroundDamp(1,:,ittm);
    else
        TBdamp_ittm.TingroundDamp = TBdamp.TingroundDamp(ittn-1,:,ittm);
    end

    
    SWRinWsun = 0; SWRinWshd = 0; 
    Tairout = SimC.TempVec.TCanyon(ittn,1,ittm);
    qairout = SimC.Humidity.CanyonSpecific(ittn,1,ittm);
    ea      = SimC.Humidity.CanyonVapourPre(ittn,1,ittm);
    G2Roof  = SimC.Gflux.G2Roof(ittn,1,ittm);
    G2WallSun   = SimC.Gflux.G2WallSun(ittn,1,ittm);
    G2WallShade = SimC.Gflux.G2WallShade(ittn,1,ittm);

%     G2Roof  = 0;
%     G2WallSun   = 0;
%     G2WallShade = 0;


[TB,fval,exitflag]=BuildingEnergyModel.fSolver_buildingint(TempVecB_ittm,MeteoData,HumidityAtm,SWRinWsun,SWRinWshd,...
    Tairout,qairout,G2Roof,G2WallSun,G2WallShade,TBdamp_ittm,...
    Gemeotry_m,PropOpticalIndoors,ParHVAC,ParCalculation,ParThermalBulidFloor);


TempVecB.Tceiling(ittn,1)	=	TB(1,1);
TempVecB.Tinwallsun(ittn,1)	=	TB(1,2);
TempVecB.Tinwallshd(ittn,1)  =	TB(1,3);
TempVecB.Tinground(ittn,1)   =	TB(1,4);
TempVecB.Tbin(ittn,1)		=	TB(1,5);
TempVecB.qbin(ittn,1)		=	TB(1,6);


[HbuildIntc,LEbuildIntc,GbuildIntc,SWRabsBc,LWRabsBc,Tdpfloor]=BuildingEnergyModel.EBSolver_BuildingIntOUTPUT(TB,MeteoData,SWRinWsun,SWRinWshd,...
    TempVecB_ittm,Tairout,qairout,G2Roof,G2WallSun,G2WallShade,TBdamp_ittm,...
    Gemeotry_m,PropOpticalIndoors,ParHVAC,ParCalculation,ParThermalBulidFloor);


%--------------------------------------------------------------------------
% Assign outputs
TBdamp.TingroundDamp(ittn,:,ittm) = Tdpfloor;

for i=1:length(HbuildIntNames)
	HbuildInt.(cell2mat(HbuildIntNames(i)))(ittn,1,ittm)=	HbuildIntc.(cell2mat(HbuildIntNames(i)));
end

for i=1:length(LEbuildIntNames)
	LEbuildInt.(cell2mat(LEbuildIntNames(i)))(ittn,1,ittm)=	LEbuildIntc.(cell2mat(LEbuildIntNames(i)));
end

for i=1:length(GbuildIntNames)
	GbuildInt.(cell2mat(GbuildIntNames(i)))(ittn,1,ittm)=	GbuildIntc.(cell2mat(GbuildIntNames(i)));
end

for i=1:length(SWRabsBNames)
	SWRabsB.(cell2mat(SWRabsBNames(i)))(ittn,1,ittm)=	SWRabsBc.(cell2mat(SWRabsBNames(i)));
end

for i=1:length(LWRabsBNames)
	LWRabsB.(cell2mat(LWRabsBNames(i)))(ittn,1,ittm)=	LWRabsBc.(cell2mat(LWRabsBNames(i)));
end


end



figure
plot(SimC.MeteoDataRaw.Date(1:n),TempVecB.Tceiling-273.15,'DisplayName','Ceiling')
hold on
plot(SimC.MeteoDataRaw.Date(1:n),TempVecB.Tinwallsun-273.15,'DisplayName','Wall sun')
plot(SimC.MeteoDataRaw.Date(1:n),TempVecB.Tinwallshd-273.15,'DisplayName','Wall shade')
plot(SimC.MeteoDataRaw.Date(1:n),TempVecB.Tinground-273.15,'DisplayName','Ground')
plot(SimC.MeteoDataRaw.Date(1:n),TempVecB.Tbin-273.15,'DisplayName','Indoor air')
plot(SimC.MeteoDataRaw.Date(1:n),SimC.TempVec.TCanyon(1:n,1,1)-273.15,'DisplayName','Outdoor air')
legend

figure
plot(SimC.MeteoDataRaw.Date(1:n),TempVecB.qbin,'DisplayName','Indoor air q')
hold on
plot(SimC.MeteoDataRaw.Date(1:n),SimC.Humidity.CanyonSpecific(1:n,1,1),'DisplayName','Outdoor air')
legend

figure
tiledlayout(2,3);
nexttile
hold on
% plot(SimC.MeteoDataRaw.Date(1:n),HbuildInt.HBinRoof,'DisplayName','Hb_{roof}')
% plot(SimC.MeteoDataRaw.Date(1:n),HbuildInt.HbinWallSun,'DisplayName','Hb_{w,sun}')
% plot(SimC.MeteoDataRaw.Date(1:n),HbuildInt.HbinWallshd,'DisplayName','Hb_{w,shade}')
% plot(SimC.MeteoDataRaw.Date(1:n),HbuildInt.HBinGround,'DisplayName','Hb_{ground}')
plot(SimC.MeteoDataRaw.Date(1:n),HbuildInt.HbuildInSurf,'DisplayName','Hb_{surftot}')
plot(SimC.MeteoDataRaw.Date(1:n),HbuildInt.Hvent,'DisplayName','Hb_{ventilation}')
plot(SimC.MeteoDataRaw.Date(1:n),HbuildInt.H_AC,'DisplayName','Hb_{AC}')
plot(SimC.MeteoDataRaw.Date(1:n),HbuildInt.H_air,'DisplayName','Hb_{air,storage}')
legend


figure
plot(SimC.MeteoDataRaw.Date,SimC.Gflux.G1WallSun(:,1,1),'DisplayName','G1_{w,sun}')
hold on
plot(SimC.MeteoDataRaw.Date,SimC.Gflux.G2WallSun(:,1,1),'DisplayName','G2_{w,sun}')
plot(SimC.MeteoDataRaw.Date,SimC.Gflux.G1WallShade(:,1,1),'DisplayName','G1_{w,shade}')
plot(SimC.MeteoDataRaw.Date,SimC.Gflux.G2WallShade(:,1,1),'DisplayName','G2_{w,shade}')
legend

Hbuild  =   10;
Wroof   =   3;

Test = Hbuild.*(SimC.Gflux.G2WallSun(:,1,1)+SimC.Gflux.G2WallShade(:,1,1))+Wroof.*SimC.Gflux.G2Roof(:,1,1);
Test2 = HbuildInt.H_AC-Test(1:n);

figure
plot(Test2)

