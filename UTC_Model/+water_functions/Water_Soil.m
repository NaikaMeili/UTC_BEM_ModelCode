function[V,O,OS,Lk,Psi_s_H,Psi_s_L,Exwat_H,Exwat_L,Rd,TE_L,TE_H,...
	E_soil,dV_dt,WBalance_soil,Psi_soil,Ko]=...
	Water_Soil(Otm1,f,TE_H,TE_L,E_soil,Qlat_in,...
	dth,Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,...
	Rrootl_H,Rrootl_L,PsiL50_H,PsiL50_L,PsiX50_H,PsiX50_L,Zs,row)

%% OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V					=	Volume of water in each soil layer [mm]. Never includes the residual water content Ohy
% O					=	Soil moisture in each soil layer [-]. Always includes the residual water content Ohy
% Lk				=	Leakage at bedrock [mm/h]
% Rd				=	saturation excess runoff / Dunne Runoff [mm] 
% Psi_s				=	Soil Water Potential for Vegetation [MPa] 
% Exwat				=	Max extractable water for Vegetation [mm m2 / m2 ground h ] 
% Psi_soil			=	Soil water potential at soil moisture O [mm]
% Ko				=	Hydraulic conductivity at soil moisture O. [mm/h}

%% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f_ground			=	Infiltration into first soil layer [mm/h]
% E_soil			=	evaporation from soil under vegetation [kg/m^2 s]
% T_soil			=	Transpiration from vegetation [kg/m^2 s]
% Otm1				=	Soil moisture in the different soil layers [-]

% Calculate soil parameters depending on soil composition
[Zs,dz,ms,Osat,Ohy,nVG,alpVG,Ks_Zs,L,Pe,O33,SPAR,EvL_Zs,Inf_Zs,...
	RfH_Zs,RfL_Zs,~,Kbot,Slo_pot,Dz,aR,aTop,~,~,~,~]...
	=soil_functions.SoilParametersTotal(Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,Zs);

E_soil			=	E_soil*3600*1000/row;	% [mm/h]
TE_L			=	TE_L*3600*1000/row;		% [mm/h]
TE_H			=	TE_H*3600*1000/row;		% [mm/h]

% Distributed Sink: Evaporation from bare soil and Transpiration
E_soil_dis		=	E_soil.*EvL_Zs;	% [mm/h] Evaporation from bare soil * ratio of soil depth
TE_dis_H		=	TE_H.*RfH_Zs;	% [mm/h] Root water uptake from different soil layers
TE_dis_L		=	TE_L.*RfL_Zs;	% [mm/h] Root water uptake from different soil layers

% Leakage bottom [mm/h]
[Lk]=soil_functions.Leakage_Bottom(Otm1,Ks_Zs,Osat,Ohy,L,nVG,Kbot,ms,SPAR); %%% [mm/h]

Lk		=	Lk;				% [mm/h]
f		=	f;				% [mm/h]
Qlat_in	=	Qlat_in/dth;	% [mm/h]

% Solving Richards equation and calculating the change in water volume
% Initial Condition --->  V- Values previous time step
% Vout is the actual water volume in the ground. The ODE solves the
% differential equation and therefore as a start we give a Volume and get
% back the volume. For the solver we just write the differential equation.
V0		=	(Otm1-Ohy).*dz;  % Water volume in each layer [mm]
T_SPAN	=	[0 dth];
ISeep	=	ones(1,ms) ;
OPT_SM	=	odeset('AbsTol',0.05,'MaxStep',dth);

[Tout,Vout]=ode23s(@soil_functions.SOIL_MOISTURES_RICH_COMP,T_SPAN,V0,OPT_SM,...
	Lk,f,E_soil_dis,TE_dis_H,TE_dis_L,Qlat_in,...
	Slo_pot,ISeep,int32(SPAR),Osat,Ohy,O33,dz,Ks_Zs,Dz,int32(ms),L,Pe,aR,aTop,alpVG,nVG,Zs,1,0,0);

% [dV]=SOIL_MOISTURES_RICH_COMP(t,V,...
% Lk,f,EG,T_H,T_L,...
% Qi_in,Slo_pot,IS,SPAR,Osat,Ohy,O33,dz,Ks_Zs,Dz,numn,L,Pe,aR,aT,alpVG,nVG,Zs,cosalp,sinalp,SN) 

V		=	Vout(end,:); % Water content in each soil layer [mm]. V never includes the residual water content Ohy

if isnan(sum(V))
	disp('NaN values in the Volumes')
	return
end

% Update soil moisture content in different soil layers.
[O,~,~,OS,Psi_s_H,Psi_s_L,~,~,Exwat_H,Exwat_L,Rd,WTR,~,~,~]=soil_functions.Soil_Water_MultiLayer(V,Zs,...
    dz,ms,Osat,Ohy,nVG,alpVG,Ks_Zs,L,Pe,O33,SPAR,EvL_Zs,Inf_Zs,RfH_Zs,RfL_Zs,Rrootl_H,Rrootl_L,PsiL50_H,PsiL50_L,PsiX50_H,PsiX50_L);
% Exwat_H, Exwat_L [mm/h]
% Rd [mm/dth]


for i=1:length(O)
[Ko(i),Psi_soil(i)]=soil_functions.Conductivity_Suction(SPAR,Ks_Zs(i),Osat(i),Ohy(i),L(i),Pe(i),O33(i),alpVG(i),nVG(i),O(i));
end

% Rd = saturation excess runoff [mm]
% V(i) = water volume in soil layer i
% WTR(i) = water flow due to water table rising
% Volume Correction for Rd and WTR
V(1)		=	V(1) +  WTR(2) - Rd;
V(2:end-1)	=	V(2:end-1)+ (WTR(3:end) - WTR(2:end-1));
V(end)		=	V(end) - WTR(end);

% Volume Compensation -- Negative Value
if sum(V < 0) > 0
	[V,TE_dis_H,TE_dis_L,E_soil,Lk]=soil_functions.Volume_Correction(V,EvL_Zs,RfH_Zs,RfL_Zs,dth*E_soil,dth*TE_dis_H,dth*TE_dis_L,dth*Lk);
	TE_H	=	sum(TE_dis_H)/dth;	% [mm/h]
	TE_L	=	sum(TE_dis_L)/dth;	% [mm/h]
	E_soil	=	E_soil/dth;			% [mm/h]
 	Lk		=	Lk/dth;				% [mm/h]
end

% Volume Balance check
dV_dt			=	sum(V,2)-sum(V0,2);
WBalance_soil	=	dth*f + dth*nansum(Qlat_in) - dth*Lk - dth*TE_H - dth*TE_L - dth*E_soil - Rd - dV_dt;	% [mm/dth]

%-----------------------------------------
E_soil			=	E_soil/(3600*1000/row);		% [kg/m^2 s]
TE_L			=	TE_L/(3600*1000/row);		% [kg/m^2 s]
TE_H			=	TE_H/(3600*1000/row);		% [kg/m^2 s]
Lk				=	Lk;							% [mm/h]
Exwat_H			=	Exwat_H;					% [mm/h]
Exwat_L			=	Exwat_L;					% [mm/h]
Rd				=	Rd;							% [mm/dth]
dV_dt			=	dV_dt;						% [mm/dth]


end



