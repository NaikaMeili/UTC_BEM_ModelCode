function[]=UrbanClimateVariables(TempVec,UTCI,Results2m,MeteoDataRaw,MeanRadiantTemperature,...
    FractionsGround_Out,FractionsRoof_Out,ParTree_Out,Figure)


TTUrban = table(hour(MeteoDataRaw.Date),month(MeteoDataRaw.Date),...
            Results2m.T2m-273.15, Results2m.RH_T2m.*100, UTCI, MeanRadiantTemperature.Tmrt,...
            TempVec.TRoofImp-273.15,TempVec.TRoofVeg-273.15,...
            TempVec.TGroundImp-273.15,TempVec.TGroundBare-273.15,TempVec.TGroundVeg-273.15,...
            TempVec.TTree-273.15,TempVec.TWallSun-273.15,TempVec.TWallShade-273.15,...
            TempVec.Tatm-273.15,MeteoDataRaw.rel_humidity.*100);

TTUrban.Properties.VariableNames = {'Hour','Month','T2m','RH2m','UTCI','Tmrt',...
    'TroofImp','TroofVeg','TgroundImp','TgroundBare','TgroundVeg','Ttree',...
    'TWallSun','TWallShade','Tatm','RHatm'};


TTUrbanDiurnal = varfun(@nanmean,TTUrban,'GroupingVariables','Hour');
TTUrbanSeasonal = varfun(@nanmean,TTUrban,'GroupingVariables','Month');


if FractionsRoof_Out.fimp>0; CFrimp = 1; else CFrimp = NaN; end
if FractionsRoof_Out.fveg>0; CFrveg = 1; else CFrveg = NaN; end
if FractionsGround_Out.fimp>0; CFgimp = 1; else CFgimp = NaN; end
if FractionsGround_Out.fbare>0; CFgbare = 1; else CFgbare = NaN; end
if FractionsGround_Out.fveg>0; CFgveg = 1; else CFgveg = NaN; end
if ParTree_Out.trees>0; CFtree = 1; else CFtree = NaN; end


if Figure==1
% Outdoor thermal comfort
%--------------------------------------------------------------------------
% Plot figures
f1 = figure;
set(f1, 'Units','centimeters','Position', [1 1 14 14])

t = tiledlayout(2,2);
t.Padding = 'compact'; %t.TileSpacing = 'compact';

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_UTCI,'r','LineWidth',1.5,'DisplayName','UTCI')
hold on
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_T2m,'b--','LineWidth',1.5,'DisplayName','T_{2m}')
xlim([0 23]); xlabel('hour'); ylabel('T (\circC)'); title('Outdoor thermal comfort'); subtitle('Diurnal');
legend('Location','southoutside','NumColumns',2)

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_Tmrt,'r','LineWidth',1.5,'DisplayName','T_{mrt}')
hold on
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_T2m,'b--','LineWidth',1.5,'DisplayName','T_{2m}')
xlim([0 23]); xlabel('hour'); ylabel('T (\circC)'); title('Mean radiant temperature'); subtitle('Diurnal');
legend('Location','southoutside','NumColumns',2)


nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_UTCI,'r','LineWidth',1.5,'DisplayName','UTCI')
hold on
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_T2m,'b--','LineWidth',1.5,'DisplayName','T_{2m}')
xlim([1 12]); xlabel('hour'); ylabel('T (\circC)'); subtitle('Seasonal');
%legend('Location','southoutside','NumColumns',2)

nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_Tmrt,'r','LineWidth',1.5,'DisplayName','T_{mrt}')
hold on
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_T2m,'b--','LineWidth',1.5,'DisplayName','T_{2m}')
xlim([1 12]); xlabel('hour'); ylabel('T (\circC)'); subtitle('Seasonal');
%legend('Location','southoutside','NumColumns',2)


% Plot figures
%--------------------------------------------------------------------------
f1 = figure;
set(f1, 'Units','centimeters','Position', [1 1 19 14])

t = tiledlayout(2,3);
t.Padding = 'compact'; %t.TileSpacing = 'compact';

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_TroofImp.*CFrimp,'r','LineWidth',1.5,'DisplayName','T_{roof,imp}')
hold on
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_TroofVeg.*CFrveg,'b','LineWidth',1.5,'DisplayName','T_{roof,veg}')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_TgroundImp.*CFgimp,'r--','LineWidth',1.5,'DisplayName','T_{ground,imp}')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_TgroundBare.*CFgbare,'y','LineWidth',1.5,'DisplayName','T_{ground,bare}')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_TgroundVeg.*CFgveg,'b--','LineWidth',1.5,'DisplayName','T_{ground,veg}')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_Ttree.*CFtree,'g','LineWidth',1.5,'DisplayName','T_{tree}')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_TWallSun,'m','LineWidth',1.5,'DisplayName','T_{wall,sun}')
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_TWallShade,'c','LineWidth',1.5,'DisplayName','T_{wall,shade}')
xlim([0 23]); xlabel('hour'); ylabel('T (\circC)'); title('Surface temperature'); subtitle('Diurnal');
legend('Location','southoutside','NumColumns',2)


nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_T2m,'b','LineWidth',1.5,'DisplayName','T_{air,2m}')
hold on
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_Tatm,'k','LineWidth',1.5,'DisplayName','T_{air,atm}')
xlim([0 23]); xlabel('hour'); ylabel('T (\circC)'); title('Air temperature'); subtitle('Diurnal');
legend('Location','southoutside','NumColumns',2)


nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_RH2m,'g','LineWidth',1.5,'DisplayName','RH_{2m}')
hold on
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_RHatm,'k','LineWidth',1.5,'DisplayName','RH_{atm}')
xlim([0 23]); xlabel('hour'); ylabel('T (\circC)'); title('Relative humidity'); subtitle('Diurnal');
legend('Location','southoutside','NumColumns',2)

%--------------------------------------------------------------------------
nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_TroofImp.*CFrimp,'r','LineWidth',1.5,'DisplayName','T_{roof,imp}')
hold on
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_TroofVeg.*CFrveg,'b','LineWidth',1.5,'DisplayName','T_{roof,veg}')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_TgroundImp.*CFgimp,'r--','LineWidth',1.5,'DisplayName','T_{ground,imp}')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_TgroundBare.*CFgbare,'y','LineWidth',1.5,'DisplayName','T_{ground,bare}')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_TgroundVeg.*CFgveg,'b--','LineWidth',1.5,'DisplayName','T_{ground,veg}')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_Ttree.*CFtree,'g','LineWidth',1.5,'DisplayName','T_{tree}')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_TWallSun,'m','LineWidth',1.5,'DisplayName','T_{wall,sun}')
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_TWallShade,'c','LineWidth',1.5,'DisplayName','T_{wall,shade}')
xlim([1 12]); xlabel('Month'); ylabel('T (\circC)'); subtitle('Seasonal');
%legend('Location','southoutside','NumColumns',2)


nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_T2m,'b','LineWidth',1.5,'DisplayName','T_{air,2m}')
hold on
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_Tatm,'k','LineWidth',1.5,'DisplayName','T_{air,atm}')
xlim([1 12]); xlabel('Month'); ylabel('T (\circC)'); subtitle('Seasonal');
%legend('Location','southoutside','NumColumns',2)


nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_RH2m,'g','LineWidth',1.5,'DisplayName','RH_{2m}')
hold on
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_RHatm,'k','LineWidth',1.5,'DisplayName','RH_{atm}')
xlim([1 12]); xlabel('Month'); ylabel('T (\circC)');subtitle('Seasonal');
%legend('Location','southoutside','NumColumns',2)
end
