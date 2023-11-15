function[]=UrbanClimateVariables(TempVec,UTCI,Results2m,MeteoDataRaw,MeanRadiantTemperature,...
    TempVecB,HumidityBuilding,...
    FractionsGround_Out,FractionsRoof_Out,ParTree_Out,Figure,ittmTot,BEM_on)

for ittm=1:ittmTot
TTUrban = table(hour(MeteoDataRaw.Date),month(MeteoDataRaw.Date),...
            Results2m.T2m(:,1,ittm)-273.15, TempVec.TCanyon(:,1,ittm)-273.15, ...
            Results2m.RH_T2m(:,1,ittm).*100, UTCI(:,1,ittm), MeanRadiantTemperature.Tmrt(:,1,ittm),...
            TempVec.TRoofImp(:,1,ittm)-273.15,TempVec.TRoofVeg(:,1,ittm)-273.15,...
            TempVec.TGroundImp(:,1,ittm)-273.15,TempVec.TGroundBare(:,1,ittm)-273.15,TempVec.TGroundVeg(:,1,ittm)-273.15,...
            TempVec.TTree(:,1,ittm)-273.15,TempVec.TWallSun(:,1,ittm)-273.15,TempVec.TWallShade(:,1,ittm)-273.15,...
            TempVec.Tatm(:,1,ittm)-273.15,MeteoDataRaw.rel_humidity.*100,...
            TempVecB.Tbin(:,1,ittm)-273.15,100.*HumidityBuilding.RHbin(:,1,ittm));

TTUrban.Properties.VariableNames = {'Hour','Month','T2m','Tcan','RH2m','UTCI','Tmrt',...
    'TroofImp','TroofVeg','TgroundImp','TgroundBare','TgroundVeg','Ttree',...
    'TWallSun','TWallShade','Tatm','RHatm','Tbin','RHbin'};


TTUrbanDiurnal = varfun(@nanmean,TTUrban,'GroupingVariables','Hour');
TTUrbanSeasonal = varfun(@nanmean,TTUrban,'GroupingVariables','Month');


if FractionsRoof_Out(ittm).fimp>0; CFrimp = 1; else CFrimp = NaN; end
if FractionsRoof_Out(ittm).fveg>0; CFrveg = 1; else CFrveg = NaN; end
if FractionsGround_Out(ittm).fimp>0; CFgimp = 1; else CFgimp = NaN; end
if FractionsGround_Out(ittm).fbare>0; CFgbare = 1; else CFgbare = NaN; end
if FractionsGround_Out(ittm).fveg>0; CFgveg = 1; else CFgveg = NaN; end
if ParTree_Out(ittm).trees>0; CFtree = 1; else CFtree = NaN; end


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
legend('Location','southoutside','NumColumns',2); grid on;

nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_Tmrt,'r','LineWidth',1.5,'DisplayName','T_{mrt}')
hold on
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_T2m,'b--','LineWidth',1.5,'DisplayName','T_{2m}')
xlim([0 23]); xlabel('hour'); ylabel('T (\circC)'); title('Mean radiant temperature'); subtitle('Diurnal');
legend('Location','southoutside','NumColumns',2); grid on;


nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_UTCI,'r','LineWidth',1.5,'DisplayName','UTCI')
hold on
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_T2m,'b--','LineWidth',1.5,'DisplayName','T_{2m}')
xlim([1 12]); xlabel('hour'); ylabel('T (\circC)'); subtitle('Seasonal'); grid on;
%legend('Location','southoutside','NumColumns',2)

nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_Tmrt,'r','LineWidth',1.5,'DisplayName','T_{mrt}')
hold on
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_T2m,'b--','LineWidth',1.5,'DisplayName','T_{2m}')
xlim([1 12]); xlabel('hour'); ylabel('T (\circC)'); subtitle('Seasonal'); grid on;
%legend('Location','southoutside','NumColumns',2)
sgtitle(['ittm = ' num2str(ittm)]) 


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
legend('Location','southoutside','NumColumns',2); grid on;


nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_T2m,'b','LineWidth',1.5,'DisplayName','T_{air,2m}')
hold on
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_Tcan,'b--','LineWidth',1.5,'DisplayName','T_{air,can}')
if BEM_on==1
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_Tbin,'c','LineWidth',1.5,'DisplayName','T_{b,indoors}')
end
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_Tatm,'k','LineWidth',1.5,'DisplayName','T_{air,atm}')
xlim([0 23]); xlabel('hour'); ylabel('T (\circC)'); title('Air temperature'); subtitle('Diurnal');
legend('Location','southoutside','NumColumns',2); grid on;


nexttile
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_RH2m,'g','LineWidth',1.5,'DisplayName','RH_{2m}')
hold on
if BEM_on==1
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_RHbin,'y','LineWidth',1.5,'DisplayName','RH_{b,indoors}')
end
plot(TTUrbanDiurnal.Hour,TTUrbanDiurnal.nanmean_RHatm,'k','LineWidth',1.5,'DisplayName','RH_{atm}')
xlim([0 23]); xlabel('hour'); ylabel('Relative humidity (%)'); title('Relative humidity'); subtitle('Diurnal');
legend('Location','southoutside','NumColumns',2); grid on;

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
xlim([1 12]); xlabel('Month'); ylabel('T (\circC)'); subtitle('Seasonal'); grid on;
%legend('Location','southoutside','NumColumns',2)


nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_T2m,'b','LineWidth',1.5,'DisplayName','T_{air,2m}')
hold on
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_Tcan,'b--','LineWidth',1.5,'DisplayName','T_{air,can}')
if BEM_on==1
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_Tbin,'c','LineWidth',1.5,'DisplayName','T_{b,indoors}')
end
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_Tatm,'k','LineWidth',1.5,'DisplayName','T_{air,atm}')
xlim([1 12]); xlabel('Month'); ylabel('T (\circC)'); subtitle('Seasonal'); grid on;
%legend('Location','southoutside','NumColumns',2)


nexttile
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_RH2m,'g','LineWidth',1.5,'DisplayName','RH_{2m}')
hold on
if BEM_on==1
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_RHbin,'y','LineWidth',1.5,'DisplayName','RH_{b,indoors}')
end
plot(TTUrbanSeasonal.Month,TTUrbanSeasonal.nanmean_RHatm,'k','LineWidth',1.5,'DisplayName','RH_{atm}')
xlim([1 12]); xlabel('Month'); ylabel('Relative humidity (%)');subtitle('Seasonal');
%legend('Location','southoutside','NumColumns',2)
sgtitle(['ittm = ' num2str(ittm)]) ; grid on;
end
end
