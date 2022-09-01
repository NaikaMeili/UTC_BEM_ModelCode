%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Leaf_Boundary_Resistence   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%REFERENCES %%   - Vesala 1998 -- Chodhury and Monteith 1988 - Ivanov 2008
%%% Leigh et al 2011 ; Leuning et al 1995
function[rb]=Leaf_Boundary_Resistence(Ws,Ts,Ta,hc,d_leaf,LAI,zatm,disp_h,zom)
%%%INPUTS
% Ta = %% air temperature [°C] -- 
% Ts =  ; %% surface temperature [°C] -- 
% Pre = %pressure [Pa]-- 
% z= %  reference height [m] % --- 
% zoh = %% roughness  eddy diffusivities for heat  [m]  
u = Ws; %% [m/s] wind speed  --- 
% hc = canpy height [m]
% d_leaf = [cm] Leaf Dimension % zom = roughness eddy diffusivities for momentum [m]

% LAI = [Leaf Area Index ]
%%% OUTPUTS 
%%% rb % [s/m] Leaf Boundary Resistence 
%%% PARAMETERS
d_leaf = d_leaf/100; %% [m] 
k= 0.4; %% Von Karman Constant 
a = 0.005; %% [m/s^0.5] --  Chodhury and Monteith 1988 
d = disp_h; %% Zero plane displacement [m] 
z= zatm ; %% Measurement Height [m] 
%%% Hypothesis Logaritmic distribution of wind speed 
us =  k*u/log((z-d)/zom); %%% Friction Velocity  [m/s]
u_hc = (us/k)*log((hc-d)/zom); %% Wind Speed top Canopy [m/s]
%%% 
alpha = log(u/u_hc)/(z/hc -1); %% Attenuation Coefficient  
alpha = 0.5*alpha*LAI/2; 
%%% u(z) = u_hc*exp(alpha*((z/hc-1)); 
%%%%%%%%% Expression of Leaf Boundary Layer Resistence  
%%% gb(z)= a*(u(z)/d_leaf)^0.5; %% [Jones 1992]  Integral between 0-LAI_TOT with LAI(z) = LAI_TOT*z/Hc   
gb = (2*a/alpha)*((u_hc/d_leaf)^0.5)*(1-exp(-alpha/2)); %% [m/s] 
%rb= (1/(gb*LAI)); %%   Leaf Boundary Layer Resistnce [s/m] Plant 
%rb=1/gb; %%   Leaf Boundary Layer Resistnce [s/m] one-sided for unit leaf 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Expression for free convection  (Leuning 1995, Monteith 1973) 
Dh= 1.9*1e-5; %% [m2/s] 
Gr= 1.6*1e+8*(Ts-Ta)*(d_leaf^3).*(Ts>Ta); %% [-] 
gb_free = 0.5*Dh*Gr^(0.25)/d_leaf; %%[m/s] 
%%%%%
gb=gb+gb_free; 
rb=1/gb; %%   Leaf Boundary Layer Resistnce [s/m] one-sided for unit leaf 
return 
