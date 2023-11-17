%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute Soil Thermal Properties     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[lanS,cv_Soil,CTt]=Soil_Thermal_properties(Tdp,rsd,lan_dry,lan_s,cv_s,Osat,Ohy,O)
%%REFERENCES 
%%% Noilhan and Mahfouf 1996 --- ISBA 2004
%%% Farouki (1981) --- Oleson et al., 2004 - Ivanov et al., 2008 --
%%% Boone et al., 2000 -- Johansen 1975

%%%INPUTS
% Tdp = dampening temperature
% rsd = density dry soil [kg/m^3]
% lan_dry = Thermal conductivity dry soil [W/m K]
% lan_s = Thermal conductivity soil solid [W/m K]
% cv_s = Volumetric heat capacity soil solid [J/m^3 K]
% Osat = Water content at saturation [-]
% Ohy = Residual / Hygroscopic / Wilting Point water content [-]
% O = Soil moisture, soil water content at previous time step [-]

%%% OUTPUTS
% CTt,  %  [K m^2/J] Total Thermal Capacity Soil
% lanS,  % [W/m K ] Thermal conductivity Soil 
% cv_Soil %  Volumetric heat capacity Soil  [J/m^3 K] 


%%%%%%%%%% THERMAL PROPERTIES SOIL
%%%% WATER PROPERTIES %%%%
row=1000; %%[kg/m^3]
lan_wat = 0.58; %%% [W/m K ] Thermal conductivity water
lan_ice = 2.29; %%%  [W/m K ] Thermal conductivity ice
cv_w =  4186000;  % [J/m^3 K] Volumetric heat capacity water

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=length(O);
lanS = zeros(1,n);  %  [W/m K ] Thermal conductivity Soil
cv_Soil=zeros(1,n);  %  Volumetric heat capacity Soil  [J/m^3 K]
rsoil=zeros(1,n);  %   Soil Density [kg/m^3]
cs_Soil=zeros(1,n);  %  %%% [J/kg K]  %% Specific Heat Soil
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
O(O<Ohy)=Ohy(O<Ohy); 
O(O>Osat)=Osat(O>Osat);
Oice = Osat.*(Tdp < 0); %% Frozen layers 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Each Soil layer
for i=1:n
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
	%    if Tdp(i) > 0, edited by Naika
    if Tdp > 0
        lan_sat = (lan_wat^Osat(i))*(lan_s(i)^(1-Osat(i))); %% Saturated Conductivity [W/m K]
        Ke = log((O(i)+Oice(i))/Osat(i))+ 1; %% Kersten number
        Ke= Ke*(Ke>=0);
    else
        Oliq= Osat(i) - Oice(i); %%% Liquid water content at saturation
        lan_sat = (lan_wat^Osat(i))*(lan_s(i)^(1-Osat(i)))*(lan_ice^(Osat(i)-Oliq)); %% Saturated Conductivity [W/m K]
        Ke = (O(i)+Oice(i))/Osat(i);
    end
    
    if O(i)/Osat(i) > 1*10^-7
        lanS(i) = Ke*lan_sat + (1-Ke)*lan_dry(i) ; % [W/m K ] Thermal conductivity Soil
    else
        lanS(i) = lan_dry(i);%   [W/m K ] Thermal conductivity Soil
    end
    
    %%%%%%%%%%%%%%%%%%%%%%
    cv_Soil(i) = cv_s(i)*(1-Osat(i)) + O(i)*cv_w;  %  Volumetric heat capacity Soil  [J/m^3 K]
    rsoil(i)= rsd(i) + (O(i)-Ohy(i))*row; %% Soil Density [kg/m^3]
    cs_Soil(i) = cv_Soil(i)/rsoil(i) ; %%% [J/kg K]  %% Specific Heat Soil
end

%%%%%%%%%%%%%%%%%%%%%%%%%
tau= 86400; %% [s] time constant
CTt = 2.*(sqrt(pi./(lanS.*cs_Soil.*rsoil.*tau))); %%  [K m^2/J] Total Thermal Capacity Soil

end





