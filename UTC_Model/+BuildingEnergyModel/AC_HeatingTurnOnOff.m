function[ParHVAC,ParHVACorig]=AC_HeatingTurnOnOff(ParHVAC,TempVecB_ittm,TempVec_ittm,Humidity_ittm,MeteoData,Gemeotry_m,BEM_on)

% Initial AC parameters
ParHVACorig = ParHVAC;

% Urban geometry [m]
Hwall   = Gemeotry_m.Height_canyon;     % Building height [m]
Wroof   = Gemeotry_m.Width_roof;        % Roof width [m]

% Temperature [K] and specific humidity [kg/kg] conditions from previous time step
Tbin_tm1            = TempVecB_ittm.Tbin;               % Building internal temperature from previous timestep
qbin_tm1            = TempVecB_ittm.qbin;               % Building internal humidity from previous timestep
Tintmasstm1         = TempVecB_ittm.Tintmass;           % Ground temperature from previous timestep   

Tsurftm1 = (Wroof*(TempVecB_ittm.Tinground + TempVecB_ittm.Tceiling) + Hwall*(TempVecB_ittm.Tinwallshd + TempVecB_ittm.Tinwallsun))/(2*Wroof + 2*Hwall);

% Heating and cooling parameters
ACon                = ParHVAC.ACon;                     % Parameter specifying if AC is on (=1) or off (=0)
Heatingon           = ParHVAC.Heatingon;                % Parameter specifying if building internal heating is on and off (on = 1, off = 0)
TsetpointCooling    = ParHVAC.TsetpointCooling;         % Cooling set-point temperature [K]
TsetpointHeating    = ParHVAC.TsetpointHeating;         % Heating set-point temperature [K]
RHsetpointCooling   = ParHVAC.RHsetpointCooling./100;   % Cooling set-point humidity [-]
RHsptHeating        = NaN;                              % Heating set-point humidity, currently not used in the model

% Calculate different humidity metrics for indoor humidity regulation
%--------------------------------------------------------------------------
esat_TspCooling  =    611*exp(17.27*(TsetpointCooling-273.16)/(237.3+(TsetpointCooling-273.16)));	% Saturation vapor pressure at set-point temperature [Pa]
ea_RHspCooling   =    RHsetpointCooling.*esat_TspCooling; % Indoor vapor pressure at building aircon set-point [Pa]
q_RHspCooling    =    0.622*ea_RHspCooling/(MeteoData.Pre-0.378*ea_RHspCooling);    % Specifc humidity of air at building aircon set-point [kg/kg]


% Control in case heating and cooling set-points are the same
if TsetpointCooling == TsetpointHeating  && ACon==1 && Heatingon==1
    TsetpointHeating = TsetpointHeating-0.01;
elseif TsetpointHeating > TsetpointCooling && ACon==1 && Heatingon==1
    disp("Cooling set point is lower than heating setpoint, please double check")
end


% Parameter initiation of cooling (sensible heat) and de-humidifcation
% (latent heat) during AC mode
if ACon==1
    AC_onCool   = 1; AC_onDehum  = 1;
else
    AC_onCool   = 0; AC_onDehum  = 0;
end

%--------------------------------------------------------------------------
% Automatic switch on and off air-conditioning and heating if the indoor
% temperature is likely to exceed / fall below the set threshold
% temperature and the user specified the use of air-conditioning and
% heating

% An intial guess is based on the outdoor temperature and indoor surface
% temperatures of the previous time-step to avoid issues in the numerical
% solver
if ACon == 1
    if TempVec_ittm.TCanyon>=TsetpointCooling && Tsurftm1>=TsetpointCooling
        AC_onCool = 1;
        AC_onDehum = 1;
        if Humidity_ittm.CanyonSpecific<q_RHspCooling
            AC_onDehum = 0;
        end
    else 
        ACon = 0; 
        AC_onCool   = 0; 
        AC_onDehum  = 0;
    end
end

if Heatingon==1
    if TempVec_ittm.TCanyon<=TsetpointHeating && Tsurftm1<=TsetpointHeating
        Heatingon = 1; 
    else 
        Heatingon = 0; 
    end
end

% Updated AC on/off parameters
if BEM_on == 1
    ParHVAC.ACon        = ACon;
    ParHVAC.AC_onCool   = AC_onCool;
    ParHVAC.AC_onDehum  = AC_onDehum;
    ParHVAC.Heatingon   = Heatingon;
    ParHVAC.TsetpointCooling  = TsetpointCooling;
    ParHVAC.TsetpointHeating  = TsetpointHeating;
    ParHVAC.q_RHspCooling = q_RHspCooling;
    ParHVAC.MasterOn    = 0;
else
    ParHVAC.ACon        = 0;
    ParHVAC.AC_onCool   = 0;
    ParHVAC.AC_onDehum  = 0;
    ParHVAC.Heatingon   = 0;
    ParHVAC.TsetpointCooling  = TsetpointCooling;
    ParHVAC.TsetpointHeating  = TsetpointHeating;
    ParHVAC.q_RHspCooling = q_RHspCooling;
    ParHVAC.MasterOn    = 0;
end




