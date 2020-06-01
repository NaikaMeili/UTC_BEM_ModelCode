%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Conductivity_Suction_O  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Ko,Po]=Conductivity_Suction(SPAR,Ks,Osat,Ohy,L,Pe,O33,alpVG,nVG,O) 
%%REFERENCES %%   Saxton and Rawls 2006  
%%%INPUTS
%%% Osat [] Saturation moisture 0 kPa 
%%% L % Slope of logaritimc tension-moisture curve 
%%% Pe % Tension at air antry (bubbling pressure) [kPa]
%%% Ks  % saturation conductivty [mm/h]
%%% O33 %% 33 kPa Moisture 
%%% O Soil Moisture -- [] 
%%%dt time step [s]
%%% OUTPUTS 
%%%%%% Ko  % hydraulic conductivty at O [mm/h]
%%% Po Tension at O [mm]
%Fe  Desorption rate [mm/h]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch SPAR
    case 1
        %%% Van-Genuchten, 1980
        Se = (O-Ohy)./(Osat-Ohy);
        mVG= 1-1./nVG;
        Po = -(1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
        Ko= Ks.*((Se).^(0.5)).*(1-(1-(Se).^(1./mVG)).^mVG).^2; %%% [mm/h]
    case 2        
        gw= 9810; %% specific weight water [N/m^3]
        B=1./L;
        A= exp(log(33)+B.*log(O33)); % Coefficient of moisture tension
        %%%%%%%%%%%%%%%%%%%%
        Ko = Ks.*(O./Osat).^(3+(2./L)); %%% [mm/h]
        %Ko = Ks*power((O/Osat),(3+(2/L))); %%% [mm/h]
        if O < O33
            Psi = A.*(O.^-B); %% [kPa]
            %Psi = A*power(O,-B); %% [kPa]
        else
            Psi =33-((O-O33).*(33-Pe)./(Osat-O33));%% [kPa]
        end
        %%%%%%%%%5
        Po =1000.*1000.*Psi./(gw); %%[mm]% Tension at O
        %%%%%%%%%%%%%%%
end
return 
