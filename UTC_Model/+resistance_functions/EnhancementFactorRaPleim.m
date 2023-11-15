function[fconv,ra_enhanced,ra, LAN]=EnhancementFactorRaPleim(ra,zom,zoh,disp_h,zatm,Ws,hPBL)


% According to Pleim et al. 2007
%-------------------------------------------------------------------------
%[ra]=Aerodynamic_Resistence(Ta,Ts,Pre,zatm,disp_h,zom,zoh,Ws,ea,es);
% Backcalculate Monin-Obhukov Length
[LAN,~]=resistance_functions.BackcalculateObhukovLength(ra,zom,zoh,disp_h,zatm,Ws);

% Fraction of convective (non-local transport) according to Pleim et al. 2007
a = 7.2; % According to Holtslag et al. 1993
k = 0.4; % von Karman constant

fconv  = (1+ k^(-2/3)/(0.1*a)*(-hPBL/LAN)^(-1/3))^-1;

ra_enhanced = ra.*(1-fconv);
