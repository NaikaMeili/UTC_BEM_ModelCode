function[TempVec_ittm,Int_ittm,ExWater_ittm,Vwater_ittm,Owater_ittm,SoilPotW_ittm,...
    CiCO2Leaf_ittm,Humidity_ittm,TempDamp_ittm,Runon_ittm,Qinlat_ittm,TempVecB_ittm]...
    =CreateVariablesOfPreviousTimeStep(ittm,ittn,TempVec,...
    Int,ExWater,Vwater,Owater,SoilPotW,CiCO2Leaf,Humidity,TempDamp,Runon,Qinlat,TempVecB,Results2m)


%--------------------------------------------------------------------------
TempVecNames	=	fieldnames(TempVec);
for i=1:size(TempVecNames,1)
	if ittn==1
		TempVec_ittm.(cell2mat(TempVecNames(i)))	=	TempVec.(cell2mat(TempVecNames(i)))(1,:,ittm); 
	else
		TempVec_ittm.(cell2mat(TempVecNames(i)))	=	TempVec.(cell2mat(TempVecNames(i)))(ittn-1,:,ittm);
	end
end

if ittn==1
    TempVec_ittm.T2m = Results2m.T2m(1,:,ittm); 
else
    TempVec_ittm.T2m = Results2m.T2m(ittn-1,:,ittm); 
end


IntNames	=	fieldnames(Int);
for i=1:size(IntNames,1)
	if ittn==1
		Int_ittm.(cell2mat(IntNames(i)))	=	Int.(cell2mat(IntNames(i)))(1,:,ittm); 
	else
		Int_ittm.(cell2mat(IntNames(i)))	=	Int.(cell2mat(IntNames(i)))(ittn-1,:,ittm); 
	end
end

ExWaterNames	=	fieldnames(ExWater);
for i=1:size(ExWaterNames,1)
	if ittn==1
		ExWater_ittm.(cell2mat(ExWaterNames(i)))	=	ExWater.(cell2mat(ExWaterNames(i)))(1,:,ittm);
	else
		ExWater_ittm.(cell2mat(ExWaterNames(i)))	=	ExWater.(cell2mat(ExWaterNames(i)))(ittn-1,:,ittm);
	end
end

VwaterNames	=	fieldnames(Vwater);
for i=1:size(VwaterNames,1)
	if ittn==1
		Vwater_ittm.(cell2mat(VwaterNames(i)))	=	Vwater.(cell2mat(VwaterNames(i)))(1,:,ittm);
	else
		Vwater_ittm.(cell2mat(VwaterNames(i)))	=	Vwater.(cell2mat(VwaterNames(i)))(ittn-1,:,ittm);
	end
end

OwaterNames	=	fieldnames(Owater);
for i=1:size(OwaterNames,1)
	if ittn==1
		Owater_ittm.(cell2mat(OwaterNames(i)))	=	Owater.(cell2mat(OwaterNames(i)))(1,:,ittm);
	else
		Owater_ittm.(cell2mat(OwaterNames(i)))	=	Owater.(cell2mat(OwaterNames(i)))(ittn-1,:,ittm);
	end
end

SoilPotWNames	=	fieldnames(SoilPotW);
for i=1:size(SoilPotWNames,1)
	if ittn==1
		SoilPotW_ittm.(cell2mat(SoilPotWNames(i)))	=	SoilPotW.(cell2mat(SoilPotWNames(i)))(1,:,ittm);
	else
		SoilPotW_ittm.(cell2mat(SoilPotWNames(i)))	=	SoilPotW.(cell2mat(SoilPotWNames(i)))(ittn-1,:,ittm);
	end
end

CiCO2LeafNames	=	fieldnames(CiCO2Leaf);
for i=1:size(CiCO2LeafNames,1)
	if ittn==1
		CiCO2Leaf_ittm.(cell2mat(CiCO2LeafNames(i)))	=	CiCO2Leaf.(cell2mat(CiCO2LeafNames(i)))(1,:,ittm); 
	else
		CiCO2Leaf_ittm.(cell2mat(CiCO2LeafNames(i)))	=	CiCO2Leaf.(cell2mat(CiCO2LeafNames(i)))(ittn-1,:,ittm); 
	end
end

HumidityNames	=	fieldnames(Humidity);
for i=1:size(HumidityNames,1)
	if ittn==1
		Humidity_ittm.(cell2mat(HumidityNames(i)))	=	Humidity.(cell2mat(HumidityNames(i)))(1,:,ittm); 
	else
		Humidity_ittm.(cell2mat(HumidityNames(i)))	=	Humidity.(cell2mat(HumidityNames(i)))(ittn-1,:,ittm); 
	end
end

if ittn==1
    Humidity_ittm.q2m = Results2m.q2m(1,:,ittm); 
else
    Humidity_ittm.q2m = Results2m.q2m(ittn-1,:,ittm); 
end


TempDampNames	=	fieldnames(TempDamp);
for i=1:size(TempDampNames,1)
	if ittn==1
		TempDamp_ittm.(cell2mat(TempDampNames(i)))	=	TempDamp.(cell2mat(TempDampNames(i)))(1,:,ittm);
	else
		TempDamp_ittm.(cell2mat(TempDampNames(i)))	=	TempDamp.(cell2mat(TempDampNames(i)))(ittn-1,:,ittm);
	end
end

RunonNames	=	fieldnames(Runon);
for i=1:size(RunonNames,1)
	if ittn==1
		Runon_ittm.(cell2mat(RunonNames(i)))	=	Runon.(cell2mat(RunonNames(i)))(1,:,ittm); 
	else
		Runon_ittm.(cell2mat(RunonNames(i)))	=	Runon.(cell2mat(RunonNames(i)))(ittn-1,:,ittm); 
	end
end

QinlatNames	=	fieldnames(Qinlat);
for i=1:size(QinlatNames,1)
	if ittn==1
		Qinlat_ittm.(cell2mat(QinlatNames(i)))	=	Qinlat.(cell2mat(QinlatNames(i)))(1,:,ittm); 
	else
		Qinlat_ittm.(cell2mat(QinlatNames(i)))	=	Qinlat.(cell2mat(QinlatNames(i)))(ittn-1,:,ittm); 
	end
end

%---------------------------------------------------------
TempVecBNames	=	fieldnames(TempVecB);
for i=1:size(TempVecBNames,1)
	if ittn==1
		TempVecB_ittm.(cell2mat(TempVecBNames(i)))	=	TempVecB.(cell2mat(TempVecBNames(i)))(1,:,ittm); 
	else
		TempVecB_ittm.(cell2mat(TempVecBNames(i)))	=	TempVecB.(cell2mat(TempVecBNames(i)))(ittn-1,:,ittm);
	end
end


% TempVecNames	=	fieldnames(TempVec);
% for i=1:size(TempVecNames,1)
% 	if ittn==1 || 2 || 3 || 4 || 5
% 		TempVec_ittm5.(cell2mat(TempVecNames(i)))(ittn)	=	TempVec.(cell2mat(TempVecNames(i)))(ittn,:,ittm); 
% 	else
% 		TempVec_ittm5.(cell2mat(TempVecNames(i)))(ittn)	=	TempVec.(cell2mat(TempVecNames(i)))(ittn-1,:,ittm);
% 	end
% end






