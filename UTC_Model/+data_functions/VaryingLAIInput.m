function[LAI_TimeSeries]=VaryingLAIInput(LAIVarying,Name)

if LAIVarying==1
    LAI_TimeSeries		=	load(fullfile('+data_functions', [Name,'.mat']));
    LAI_TimeSeries.LAI_R = LAI_TimeSeries.LAI_grass;
    LAI_TimeSeries.LAI_G = LAI_TimeSeries.LAI_grass;
    LAI_TimeSeries.LAI_T = LAI_TimeSeries.LAI_dec;
else
    LAI_TimeSeries = NaN;
end