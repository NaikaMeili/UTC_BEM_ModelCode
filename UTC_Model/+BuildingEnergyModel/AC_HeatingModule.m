function[AC_on,AC_onCool,AC_onDehum,Heat_on,H_AC_Heat,LE_AC_Heat]=...
    AC_HeatingModule(AC_on,Heat_on,AC_onCool,AC_onDehum,ParHVAC,...
    HbuildIn,Hvent,Hequip,Hpeople,dSH_air,LEvent,LEequip,LEpeople,dSLE_air)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AC & heating module: Sensible and latent heat is removed or added as
% needed to keep air temperature and humidity constant for AC and air
% temperature constant for heating (humidity is not controlled during heating)
%--------------------------------------------------------------------------
% Cooling --------------------
if AC_on == 1
    if AC_onCool==1 && AC_onDehum==1
        H_AC_Heat    = HbuildIn + Hvent + Hequip + Hpeople - dSH_air;
        LE_AC_Heat   = LEvent + LEequip + LEpeople - dSLE_air;

%         if ParHVAC.MasterOn==0
%             if H_AC_Heat<0 % Prevent negative energy usage
%                 AC_onCool = 0;
%                 H_AC_Heat = 0;
%             end
%             if AC_onDehum==1 && LE_AC_Heat<0
%                 AC_onDehum = 0;
%                 LE_AC_Heat = 0;
%             end
%         end
    elseif AC_onCool==1 && AC_onDehum==0
        H_AC_Heat    = HbuildIn + Hvent + Hequip + Hpeople - dSH_air;
        LE_AC_Heat   = 0;

%         if ParHVAC.MasterOn==0
%             if H_AC_Heat<0 % Prevent negative energy usage
%                 AC_onCool = 0;
%                 H_AC_Heat = 0;
%             end
%         end

    elseif AC_onCool==0 && AC_onDehum==1
        H_AC_Heat    = 0;
        LE_AC_Heat   = LEvent + LEequip + LEpeople - dSLE_air;

    else
        H_AC_Heat    = 0;
        LE_AC_Heat   = 0;
    end
    
% Heating -------------------
elseif Heat_on==1
    H_AC_Heat    = HbuildIn + Hvent + Hequip + Hpeople - dSH_air;
    LE_AC_Heat   = 0;

%     if ParHVAC.MasterOn==0
%         if H_AC_Heat>0
%             Heat_on = 0;
%             H_AC_Heat = 0;
%         end
%     end

% Neither -------------------------
else
    H_AC_Heat    = 0;
    LE_AC_Heat   = 0;
end

end







