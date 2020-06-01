%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Leakage Bottom             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Lk]=Leakage_Bottom(O,Ks_Zs,Osat,Ohy,L,nVG,Kbot,ms,SPAR)
%%REFERENCES %%
%%%INPUTS
%%% OUTPUTS
%%% Lk %% [mm/h]
%%%%%%% Kbot [mm/h] hydraulic conductivty bedrock
if isnan(Kbot)
    switch SPAR
        case 1
            Se = (O(ms)-Ohy(ms))/(Osat(ms)-Ohy(ms)); mVG= 1-1/nVG(ms);
            Ko= Ks_Zs(ms)*((Se)^(0.5))*(1-(1-(Se)^(1/mVG))^mVG)^2; %%% [mm/h]
        case 2
            Ko = Ks_Zs(ms)*(O(ms)./Osat(ms))^(3+(2/L(ms))); %%% [mm/h] %%% Estimation Conductivty last layer
    end  
    Lk=Ko;
else
    %Lk = exp((log(Kbot) + log(Ko))/2);% [mm/h] leakage between layer n and bedrock
    if O(ms)> Osat(ms)-1e-5
        Lk = Kbot;
    else
        Lk = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end