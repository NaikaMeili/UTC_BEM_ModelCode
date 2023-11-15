function[TempVec_ittm2Ext,Humidity_ittm2Ext,TempVecB_ittm2Ext,Meteo_ittm]...
    =CreateVariablesOfPrevious2TimeStep(ittm,ittn,TempVec,Humidity,TempVecB,MeteoDataRaw)

%--------------------------------------------------------------------------
SWRin = MeteoDataRaw.SAB1_in + MeteoDataRaw.SAB1_in + MeteoDataRaw.SAD1_in + MeteoDataRaw.SAD2_in;

if ittn==1
	Meteo_ittm.SWRin	=	[SWRin(1,:,ittm); SWRin(ittn,:,ittm)]; 
    Meteo_ittm.Rain	    =	[MeteoDataRaw.rain(1,:,ittm); MeteoDataRaw.rain(ittn,:,ittm)];
else
	Meteo_ittm.SWRin	=	[SWRin(ittn-1,:,ittm); SWRin(ittn,:,ittm)]; 
    Meteo_ittm.Rain	    =	[MeteoDataRaw.rain(ittn-1,:,ittm); MeteoDataRaw.rain(ittn,:,ittm)];
end


%--------------------------------------------------------------------------
TempVecNames	=	fieldnames(TempVec);
for i=1:size(TempVecNames,1)
	if ittn==1 || ittn==2
		TempVec_ittm2.(cell2mat(TempVecNames(i)))	=	[TempVec.(cell2mat(TempVecNames(i)))(1,:,ittm); TempVec.(cell2mat(TempVecNames(i)))(1,:,ittm)]; 
	else
		TempVec_ittm2.(cell2mat(TempVecNames(i)))	=	TempVec.(cell2mat(TempVecNames(i)))(ittn-2:ittn-1,:,ittm);
	end
end



HumidityNames	=	fieldnames(Humidity);
for i=1:size(HumidityNames,1)
	if ittn==1 || ittn==2
		Humidity_ittm2.(cell2mat(HumidityNames(i)))	=	[Humidity.(cell2mat(HumidityNames(i)))(1,:,ittm); Humidity.(cell2mat(HumidityNames(i)))(1,:,ittm)]; 
	else
		Humidity_ittm2.(cell2mat(HumidityNames(i)))	=	Humidity.(cell2mat(HumidityNames(i)))(ittn-2:ittn-1,:,ittm); 
	end
end


%---------------------------------------------------------
TempVecBNames	=	fieldnames(TempVecB);
for i=1:size(TempVecBNames,1)
	if ittn==1 || ittn==2
		TempVecB_ittm2.(cell2mat(TempVecBNames(i)))	=	[TempVecB.(cell2mat(TempVecBNames(i)))(1,:,ittm); TempVecB.(cell2mat(TempVecBNames(i)))(1,:,ittm)]; 
	else
		TempVecB_ittm2.(cell2mat(TempVecBNames(i)))	=	TempVecB.(cell2mat(TempVecBNames(i)))(ittn-2:ittn-1,:,ittm);
	end
end



% Extrapolation of variables with linear regression
for i=1:size(TempVecNames,1)
    y0 = TempVec_ittm2.(cell2mat(TempVecNames(i)))';
    x0 = [1,2]; x02 = [3 4];
    X1 = [ones(length(x0),1)  x0']; b = X1\y0';
    TempVec_ittm2Ext.(cell2mat(TempVecNames(i))) = TempVec_ittm2.(cell2mat(TempVecNames(i)));
    TempVec_ittm2Ext.(cell2mat(TempVecNames(i)))(3:4) = b(1) + x02*b(2);
end

for i=1:size(HumidityNames,1)
    y0 = Humidity_ittm2.(cell2mat(HumidityNames(i)))';
    x0 = [1,2]; x02 = [3 4];
    X1 = [ones(length(x0),1)  x0']; b = X1\y0';
    Humidity_ittm2Ext.(cell2mat(HumidityNames(i))) = Humidity_ittm2.(cell2mat(HumidityNames(i)));
    Humidity_ittm2Ext.(cell2mat(HumidityNames(i)))(3:4) = b(1) + x02*b(2);
end

for i=1:size(TempVecBNames,1)
    y0 = TempVecB_ittm2.(cell2mat(TempVecBNames(i)))';
    x0 = [1,2]; x02 = [3 4]; 
    X1 = [ones(length(x0),1)  x0']; b = X1\y0';
    TempVecB_ittm2Ext.(cell2mat(TempVecBNames(i))) = TempVecB_ittm2.(cell2mat(TempVecBNames(i)));
    TempVecB_ittm2Ext.(cell2mat(TempVecBNames(i)))(3:4) = b(1) + x02*b(2);
end


% figure
% plot(x0,y0,'o')
% hold on
% plot(x02,y,'--r')

